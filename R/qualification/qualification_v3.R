# =============================================================================
# qualification_v3.R
# Phase 3.5 -- Model Qualification & Label Verification (v3)
#
# HANDOFF_TO_NARI_v3.md 7-3 "새 Table 2 (D10)" 구현.
#
# v2 대비 바뀐 것 (D10, 확정)
#   PASS/FAIL/CONCERN 판정 열을 전부 없앰. "관찰값, 모델값, 모델/관찰 비,
#   출처, 비고" 열만 쓰고, 세 부분으로 나눔:
#     (a) 보정 목표  -- 맞는 게 당연하므로 검증이 아니라고 각주에 씀
#     (b) 독립 점검  -- 보정에 쓰지 않은 라벨/문헌값과 대조
#     (c) 반감기     -- 모델 t½ 대 라벨 범위 (판정 없이 숫자만)
#   Hiemke 2018(AGNP) 치료적 참고범위 비교와 Michelson 2000 대조(96h/120h
#   잔존비율)는 환자별(N=10000) 분포가 필요해 tables_v3.R로 옮김(이
#   스크립트는 전형환자 1명, IIV 없음 -- calibration/typical_v3.py와 같은
#   종류의 결정론적 점검).
#
# Input : data/pk_parameters/pk_parameters_v3.R
# Output: outputs/qualification_v3/model_qualification_table2_v3.csv
#
# 실행: source("R/qualification/qualification_v3.R")
# =============================================================================

suppressPackageStartupMessages({
  library(rxode2)
  library(tidyverse)
})

BASE    <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR <- file.path(BASE, "outputs", "qualification_v3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "data", "pk_parameters", "pk_parameters_v3.R"))

# ── 1. ODE Models (v2와 동일 구조) ────────────────────────────────────────────

ode_1cmt <- rxode2({
  d/dt(depot) = -Ka * depot
  d/dt(Cp)    =  Ka * depot / Vd - (CL / Vd) * Cp
})

ode_parent_met <- rxode2({
  d/dt(depot) = -Ka  * depot
  d/dt(Cp)    =  Ka  * depot / Vd_p - (CL_p / Vd_p) * Cp
  d/dt(Cm)    =  Fm  * CL_p  * Cp   / Vd_m - (CL_m  / Vd_m) * Cm
  Ctotal      =  Cp + Cm
})

make_events <- function(dose_mg, tau_h, n_doses, wash_h, dt = 0.25) {
  ss_end <- n_doses * tau_h
  et(amt  = dose_mg,     # F=1 (D3) -- 전체 공칭 용량이 들어감
     ii   = tau_h,
     addl = n_doses - 1,
     time = 0) %>%
    et(seq(0, ss_end + wash_h, by = dt))
}

# 마지막 투여 간격(steady state) 안에서의 Cmax/Cmin/AUC0-tau
ss_window <- function(df, tau_h, n_doses, conc_col) {
  t0 <- (n_doses - 1) * tau_h
  t1 <- n_doses * tau_h
  df %>% filter(time >= t0, time <= t1) %>%
    transmute(t_rel = time - t0, conc = .data[[conc_col]])
}

get_cmax <- function(ss) max(ss$conc, na.rm = TRUE)
get_cmin <- function(ss) min(ss$conc[ss$t_rel > 0.01], na.rm = TRUE)  # 투여 직전(다음 투여 전) 최저값
trapz    <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
get_auc_tau <- function(ss) trapz(ss$t_rel, ss$conc)

get_thalf_sim <- function(df, conc_col = "Cp", frac = 0.35, floor = 0.005) {
  cmax <- max(df[[conc_col]], na.rm = TRUE)
  df_t <- df %>% filter(.data[[conc_col]] > cmax * floor, .data[[conc_col]] > 0)
  n <- nrow(df_t)
  if (n < 10) return(NA_real_)
  terminal <- df_t[ceiling(n * (1 - frac)):n, ]
  fit   <- lm(log(terminal[[conc_col]]) ~ terminal$time)
  slope <- coef(fit)[2]
  if (is.na(slope) || slope >= 0) return(NA_real_)
  -log(2) / slope
}

run_1cmt <- function(Ka, Vd, CL, dose_mg, tau_h, n_doses, wash_h, dt = 0.25) {
  ev <- make_events(dose_mg, tau_h, n_doses, wash_h, dt)
  rxSolve(ode_1cmt, params = c(Ka = Ka, Vd = Vd, CL = CL),
          events = ev, inits = c(depot = 0, Cp = 0)) %>% as.data.frame()
}

run_parent_met <- function(Ka, Vd_p, CL_p, Fm, Vd_m, CL_m,
                            dose_mg, tau_h, n_doses, wash_h, dt = 0.25) {
  ev <- make_events(dose_mg, tau_h, n_doses, wash_h, dt)
  rxSolve(ode_parent_met,
          params = c(Ka = Ka, Vd_p = Vd_p, CL_p = CL_p, Fm = Fm, Vd_m = Vd_m, CL_m = CL_m),
          events = ev, inits = c(depot = 0, Cp = 0, Cm = 0)) %>% as.data.frame()
}

row_ <- function(item, observed, model, source, note) {
  tibble(Item = item,
         Observed_ngmL_or_h = observed,
         Model_ngmL_or_h    = round(model, 2),
         Model_over_Obs     = round(model / observed, 3),
         Source             = source,
         Note               = note)
}

rows <- list()

# ── 2. (a) 보정 목표 -- "맞는 게 당연하므로 검증이 아님" ─────────────────────

p <- pk_params$sertraline
sim <- run_1cmt(p$Ka, p$Vd, p$CL, dose_mg = 200, tau_h = 24, n_doses = 30, wash_h = 48, dt = 0.5)
ss  <- ss_window(sim, 24, 30, "Cp")
rows$a1 <- row_("Sertraline 200 mg/day x30d, AUC0-24 (ng.h/mL)",
                2570, get_auc_tau(ss) * 1000,
                "Zoloft 2009 FDA label",
                "보정 목표(CL/F가 이 값에서 역산됨) -- 검증 아님(각주)")

p <- pk_params$paroxetine
sim <- run_1cmt(p$Ka, p$Vd, p$CL, dose_mg = 30, tau_h = 24, n_doses = 30, wash_h = 48, dt = 0.25)
ss  <- ss_window(sim, 24, 30, "Cp")
rows$a2 <- row_("Paroxetine 30 mg/day x30d, Cmax (ng/mL)",
                61.7, get_cmax(ss) * 1000,
                "Paxil 2017 FDA label",
                "보정 목표 -- 검증 아님(각주)")
rows$a3 <- row_("Paroxetine 30 mg/day x30d, Cmin (ng/mL)",
                30.7, get_cmin(ss) * 1000,
                "Paxil 2017 FDA label",
                "보정 목표 -- 검증 아님(각주)")

p <- pk_params$venlafaxine_XR; odv <- p$ODV
sim <- run_parent_met(p$Ka, p$Vd, p$CL, odv$Fm, odv$Vd_odv, odv$CL_odv,
                       dose_mg = 150, tau_h = 24, n_doses = 14, wash_h = 48, dt = 0.25)
ss  <- ss_window(sim, 24, 14, "Cp"); ss_m <- ss_window(sim, 24, 14, "Cm")
rows$a4 <- row_("Venlafaxine XR 150 mg QD, Cmax parent (ng/mL)",
                150, get_cmax(ss) * 1000,
                "Effexor XR 2017 FDA label Table 16",
                "보정 목표(D7, brentq) -- 검증 아님(각주)")
rows$a5 <- row_("Venlafaxine XR 150 mg QD, Cmax ODV (ng/mL)",
                260, get_cmax(ss_m) * 1000,
                "Effexor XR 2017 FDA label Table 16",
                "보정 목표(D7, brentq) -- 검증 아님(각주)")

# ── 3. (b) 독립 점검 -- 보정에 쓰지 않은 값과 대조 ───────────────────────────

p <- pk_params$sertraline
sim <- run_1cmt(p$Ka, p$Vd, p$CL, dose_mg = 200, tau_h = 24, n_doses = 30, wash_h = 48, dt = 0.5)
ss  <- ss_window(sim, 24, 30, "Cp")
rows$b1 <- row_("Sertraline 200 mg/day x30d, Cmax (ng/mL)",
                142, get_cmax(ss) * 1000,
                "Zoloft 2009 FDA label",
                "독립 점검(보정에 쓰지 않음 -- 보정은 AUC만 씀)")

p <- pk_params$fluoxetine; nfx <- p$norfluoxetine
sim <- run_parent_met(p$Ka, p$Vd, p$CL, nfx$Fm, nfx$Vd_met, nfx$CL_met,
                       dose_mg = 40, tau_h = 24, n_doses = 30, wash_h = 48, dt = 1)
ss_p <- ss_window(sim, 24, 30, "Cp"); ss_m <- ss_window(sim, 24, 30, "Cm")
flx_trough <- get_cmin(ss_p) * 1000
nfx_trough <- get_cmin(ss_m) * 1000
rows$b2 <- row_("Fluoxetine 40 mg/day x30d, trough (ng/mL)",
                NA_real_, flx_trough,
                "Prozac 2008 FDA label (범위 91-302)",
                "라벨 채혈시점 불명 -- 모델 최저농도(trough)와 비교(각주)")
rows$b3 <- row_("Norfluoxetine 40 mg/day x30d, trough (ng/mL)",
                NA_real_, nfx_trough,
                "Prozac 2008 FDA label (범위 72-258)",
                "라벨 채혈시점 불명 -- 모델 최저농도(trough)와 비교(각주)")
rows$b4 <- row_("NFX:FLX ratio (trough 기준)",
                0.86, nfx_trough / flx_trough,
                "Sagahon-Azua 2021 Table 1 (중앙값 0.86, IQR 0.57-1.23)",
                "환자별 비의 중앙값과 비교. 원 연구 채혈은 마지막 복용 후 중앙값 11.0h(IQR 3.9-24.7h)로 불균일 -- 모델은 trough(tau=24h) 기준이라 직접 비교엔 한계")

p <- pk_params$venlafaxine_IR; odv <- p$ODV
sim <- run_parent_met(p$Ka, p$Vd, p$CL, odv$Fm, odv$Vd_odv, odv$CL_odv,
                       dose_mg = 75, tau_h = 12, n_doses = 30, wash_h = 24, dt = 0.1)
ss_p <- ss_window(sim, 12, 30, "Cp"); ss_m <- ss_window(sim, 12, 30, "Cm")
ven_cmax <- get_cmax(ss_p) * 1000; odv_cmax <- get_cmax(ss_m) * 1000
rows$b5 <- row_("Venlafaxine IR 75 mg BID, Cmax parent (ng/mL)",
                225, ven_cmax,
                "Effexor XR 2017 FDA label Table 16",
                "독립 점검(D7 보정은 XR만 씀, IR은 쓰지 않음)")
rows$b6 <- row_("Venlafaxine IR 75 mg BID, Cmax ODV (ng/mL)",
                290, odv_cmax,
                "Effexor XR 2017 FDA label Table 16",
                "독립 점검(D7 보정은 XR만 씀, IR은 쓰지 않음)")

ven_avg <- get_auc_tau(ss_p) / max(ss_p$t_rel)
odv_avg <- get_auc_tau(ss_m) / max(ss_m$t_rel)
rows$b7 <- row_("ODV:venlafaxine ratio (평균농도, IR 75mg BID 기준)",
                NA_real_, odv_avg / ven_avg,
                "Klamerus 1992 초록 (AUC 비 2-3배)",
                "venlafaxine + ODV 평균농도 비. Klamerus는 AUC 비를 보고 -- 정상상태에서는 평균농도 비와 동일")

p <- pk_params$paroxetine
sim_nm <- run_1cmt(p$Ka, p$Vd, p$CL, dose_mg = 30, tau_h = 24, n_doses = 30, wash_h = 48, dt = 0.25)
ss_nm  <- ss_window(sim_nm, 24, 30, "Cp")
sim_pm <- run_1cmt(p$Ka, p$Vd, p$CL * cyp2d6$poor_metabolizer$CL_factor,
                    dose_mg = 30, tau_h = 24, n_doses = 30, wash_h = 48, dt = 0.25)
ss_pm  <- ss_window(sim_pm, 24, 30, "Cp")
rows$b8 <- row_("Paroxetine 30 mg/day, NM(EM) 정상상태 최저농도 (ng/mL)",
                27, get_cmin(ss_nm) * 1000,
                "Sindrup 1992 51:278-87 초록 (중앙값 81 nmol/L x0.33)",
                "독립 점검(다른 연구, 다른 대상군 -- Paxil 라벨과 별개)")
rows$b9 <- row_("Paroxetine 30 mg/day, PM 정상상태 최저농도 (ng/mL)",
                50, get_cmin(ss_pm) * 1000,
                "Sindrup 1992 51:278-87 초록 (중앙값 151 nmol/L x0.33)",
                "독립 점검. CL x0.58(D5) 적용")

sim_227 <- run_1cmt(p$Ka, p$Vd, p$CL, dose_mg = 22.7, tau_h = 24, n_doses = 30, wash_h = 48, dt = 0.25)
ss_227  <- ss_window(sim_227, 24, 30, "Cp")
rows$b10 <- row_("Paroxetine 22.7 mg/day(평균), NM 아침 투여 직전 농도 (ng/mL)",
                 34.4, get_cmin(ss_227) * 1000,
                 "Liao 2024 Table 1 (기하평균 34.4, 95%CI 31.3-37.7)",
                 "독립 점검. 모델을 Liao의 평균 용량(22.7mg/day)으로 재실행")

# ── 4. (c) 반감기: 모델 t½ 대 라벨 (판정 없이 숫자만, H7) ───────────────────

half_life_rows <- list()

p <- pk_params$sertraline
sim <- run_1cmt(p$Ka, p$Vd, p$CL, dose_mg = 50, tau_h = 24, n_doses = 14, wash_h = 300, dt = 1)
t_last <- 13 * 24
wash <- sim %>% filter(time >= t_last) %>% mutate(time = time - t_last)
half_life_rows$sert <- tibble(Drug = "Sertraline", thalf_sim_h = round(get_thalf_sim(wash, "Cp"), 1),
                               thalf_analytical_h = round(log(2) * p$Vd / p$CL, 1),
                               thalf_label = "22-36 h (FDA label 범위)", Note = "")

p <- pk_params$paroxetine
sim <- run_1cmt(p$Ka, p$Vd, p$CL, dose_mg = 30, tau_h = 24, n_doses = 14, wash_h = 300, dt = 0.5)
t_last <- 13 * 24
wash <- sim %>% filter(time >= t_last) %>% mutate(time = time - t_last)
half_life_rows$par <- tibble(Drug = "Paroxetine (NM)", thalf_sim_h = round(get_thalf_sim(wash, "Cp"), 1),
                              thalf_analytical_h = round(log(2) * p$Vd / p$CL, 1),
                              thalf_label = "7-65 h (FDA label 범위, 넓음 -- 비선형)", Note = "")

p <- pk_params$fluoxetine; nfx <- p$norfluoxetine
sim <- run_parent_met(p$Ka, p$Vd, p$CL, nfx$Fm, nfx$Vd_met, nfx$CL_met,
                       dose_mg = 20, tau_h = 24, n_doses = 60, wash_h = 2000, dt = 2)
t_last <- 59 * 24
wash <- sim %>% filter(time >= t_last) %>% mutate(time = time - t_last)
half_life_rows$flx <- tibble(Drug = "Fluoxetine", thalf_sim_h = round(get_thalf_sim(wash, "Cp"), 1),
                              thalf_analytical_h = round(log(2) * p$Vd / p$CL, 1),
                              thalf_label = "96-144 h (Prozac 만성투여 4-6일)",
                              Note = "모델 t1/2가 라벨 범위보다 깁니다 -- 그대로 보고(H7, 판정 아님)")
half_life_rows$nfx <- tibble(Drug = "Norfluoxetine", thalf_sim_h = round(get_thalf_sim(wash, "Cm"), 1),
                              thalf_analytical_h = round(log(2) * nfx$Vd_met / nfx$CL_met, 1),
                              thalf_label = "96-384 h (Prozac 만성투여 4-16일)",
                              Note = "v2 Table 2의 120-192h/192-384h 범위는 라벨 원문과 다름 -- 쓰지 않음(7-3)")

p <- pk_params$venlafaxine_XR; odv <- p$ODV
sim <- run_parent_met(p$Ka, p$Vd, p$CL, odv$Fm, odv$Vd_odv, odv$CL_odv,
                       dose_mg = 75, tau_h = 24, n_doses = 7, wash_h = 80, dt = 0.25)
t_last <- 6 * 24
wash <- sim %>% filter(time >= t_last) %>% mutate(time = time - t_last)
half_life_rows$ven <- tibble(Drug = "Venlafaxine", thalf_sim_h = round(get_thalf_sim(wash, "Cp"), 1),
                              thalf_analytical_h = round(log(2) * p$Vd / p$CL, 1),
                              thalf_label = "3-7 h (라벨, 5+-2h)", Note = "")
half_life_rows$odv <- tibble(Drug = "ODV", thalf_sim_h = round(get_thalf_sim(wash, "Cm"), 1),
                              thalf_analytical_h = round(log(2) * odv$Vd_odv / odv$CL_odv, 1),
                              thalf_label = "9-13 h (라벨, 11+-2h)",
                              Note = "모델 13.1h는 상한 13h에 근접/약간 초과 -- D7 근거 요약(4절) 참고, 판정 아님")

thalf_table <- bind_rows(half_life_rows)

# ── 5. 저장 ───────────────────────────────────────────────────────────────────

table2_ab <- bind_rows(rows) %>%
  mutate(Part = c(rep("(a) 보정 목표", 5), rep("(b) 독립 점검", 10)), .before = 1)

csv_path_ab <- file.path(OUT_DIR, "model_qualification_table2_v3.csv")
write_csv(table2_ab, csv_path_ab)

csv_path_c <- file.path(OUT_DIR, "model_qualification_table2c_halflife_v3.csv")
write_csv(thalf_table, csv_path_c)

cat("\n=== Table 2 (D10): calibration targets + independent checks (v3) ===\n")
print(table2_ab, n = Inf, width = Inf)
cat(sprintf("\nSaved: %s\n", csv_path_ab))

cat("\n=== Table 2(c): half-life, model vs label (no pass/fail, H7) ===\n")
print(thalf_table, n = Inf, width = Inf)
cat(sprintf("Saved: %s\n", csv_path_c))

cat("\nNote: Hiemke 2018 참고범위 비교와 Michelson 2000 대조(96h/120h 잔존비율)는\n")
cat("환자별(N=10000) 분포가 필요해 tables_v3.R에서 mc_results_v3_all을 읽어 만듦.\n\n")

cat("다음 단계: mc_params_v3.R -> mc_run_v3.R 본 실행(N=10000) 완료 후 tables_v3.R\n\n")
