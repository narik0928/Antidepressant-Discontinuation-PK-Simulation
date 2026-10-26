# =============================================================================
# qualification_v2.R
# Phase 3.5 — Model Qualification & Label Verification (v2)
#
# 원본 R/qualification/model_qualification.R은 건드리지 않음. 새 파일로 분리.
#
# v1 대비 바뀐 것
#   D3   : pk_parameters_v2.R을 source() 함 → F=1.0 통일 반영.
#   3-3  : **핵심 수정** — parent+metabolite 약물(Fluoxetine, Venlafaxine
#          IR/XR)에서 기존 코드는 Css_sim_ng 한 컬럼에 Ctotal(=Cp+Cm, 합)을
#          담아 "Fluoxetine"/"Venlafaxine" 행에 그대로 보고했음. 원고
#          Table 2에는 이 값이 "Fluoxetine (parent) 111.9 ng/mL"로 잘못
#          표기되어 NFX:FLX 비율이 56.6/111.9=0.51로 계산됐었음
#          (HANDOFF_TO_NARI.md 3-3; GPT_STAGE2_COMPARISON.md G6가 동일
#          숫자로 재현: parent 55.355, norfluoxetine 56.563, sum 111.918,
#          true ratio = 1.022). v2는 Cp(parent)/Cm(metabolite)/Ctotal(합)을
#          항상 세 컬럼으로 분리해서 출력 — 더 이상 하나로 합쳐서 내보내지
#          않음. parent:metabolite 비율도 직접 계산해서 같이 저장.
#
# Input : data/pk_parameters/pk_parameters_v2.R
# Output: outputs/qualification_v2/model_qualification_table_v2.csv
#         outputs/qualification_v2/model_qualification_plot_v2.png
#
# 실행: source("R/qualification/qualification_v2.R")
# =============================================================================

suppressPackageStartupMessages({
  library(rxode2)
  library(tidyverse)
})

# ── 0. Paths & seed ───────────────────────────────────────────────────────────
BASE    <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR <- file.path(BASE, "outputs", "qualification_v2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(2025)

source(file.path(BASE, "data", "pk_parameters", "pk_parameters_v2.R"))

# ── 1. Parameters (v2: F는 전부 1.0, pk_parameters_v2.R에서 읽음) ─────────────

pk <- list(

  sertraline = list(
    Ka      = pk_params$sertraline$Ka,
    Vd      = pk_params$sertraline$Vd,
    CL      = pk_params$sertraline$CL,
    F       = pk_params$sertraline$F,     # v2: 1.0 (변경 없음)
    dose_mg = 50,
    tau_h   = 24,
    n_doses = 14,
    wash_h  = 200
  ),

  fluoxetine = list(
    Ka      = pk_params$fluoxetine$Ka,
    Vd_p    = pk_params$fluoxetine$Vd,
    CL_p    = pk_params$fluoxetine$CL,
    F       = pk_params$fluoxetine$F,                        # v2: 0.70 → 1.0
    Fm      = pk_params$fluoxetine$norfluoxetine$Fm,
    Vd_m    = pk_params$fluoxetine$norfluoxetine$Vd_met,
    CL_m    = pk_params$fluoxetine$norfluoxetine$CL_met,
    dose_mg = 20,
    tau_h   = 24,
    n_doses = 60,
    wash_h  = 1500
  ),

  paroxetine = list(
    Ka      = pk_params$paroxetine$Ka,
    Vd      = pk_params$paroxetine$Vd,
    CL      = pk_params$paroxetine$CL,
    F       = pk_params$paroxetine$F,     # v2: 0.50 → 1.0
    dose_mg = 20,
    tau_h   = 24,
    n_doses = 14,
    wash_h  = 150
  ),

  venlafaxine_IR = list(
    Ka      = pk_params$venlafaxine_IR$Ka,
    Vd_p    = pk_params$venlafaxine_IR$Vd,
    CL_p    = pk_params$venlafaxine_IR$CL,
    F       = pk_params$venlafaxine_IR$F,            # v2: 0.45 → 1.0
    Fm      = pk_params$venlafaxine_IR$ODV$Fm,
    Vd_m    = pk_params$venlafaxine_IR$ODV$Vd_odv,
    CL_m    = pk_params$venlafaxine_IR$ODV$CL_odv,
    dose_mg = 25,
    tau_h   = 8,
    n_doses = 21,
    wash_h  = 80
  ),

  venlafaxine_XR = list(
    Ka      = pk_params$venlafaxine_XR$Ka,
    Vd_p    = pk_params$venlafaxine_XR$Vd,
    CL_p    = pk_params$venlafaxine_XR$CL,
    F       = pk_params$venlafaxine_XR$F,            # v2: 0.45 → 1.0
    Fm      = pk_params$venlafaxine_XR$ODV$Fm,
    Vd_m    = pk_params$venlafaxine_XR$ODV$Vd_odv,
    CL_m    = pk_params$venlafaxine_XR$ODV$CL_odv,
    dose_mg = 75,
    tau_h   = 24,
    n_doses = 7,
    wash_h  = 80
  )
)

# ── 2. ODE Model Definitions (v1과 동일) ──────────────────────────────────────

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

# ── 3. Helper Functions (v1과 동일) ───────────────────────────────────────────

get_tmax_cmax <- function(df, conc_col = "Cp") {
  idx  <- which.max(df[[conc_col]])
  list(Tmax_h = round(df$time[idx], 2),
       Cmax   = round(df[[conc_col]][idx], 5))
}

get_thalf_sim <- function(df, conc_col = "Cp",
                          frac = 0.35, floor = 0.005) {
  cmax <- max(df[[conc_col]], na.rm = TRUE)
  df_t <- df %>% filter(.data[[conc_col]] > cmax * floor,
                        .data[[conc_col]] > 0)
  n <- nrow(df_t)
  if (n < 10) return(NA_real_)
  terminal  <- df_t[ceiling(n * (1 - frac)) : n, ]
  fit       <- lm(log(terminal[[conc_col]]) ~ terminal$time)
  slope     <- coef(fit)[2]
  if (is.na(slope) || slope >= 0) return(NA_real_)
  round(-log(2) / slope, 1)
}

get_css <- function(df, conc_col = "Cp") {
  mean(df[[conc_col]], na.rm = TRUE)
}

make_events <- function(dose_F, tau_h, n_doses, wash_h, dt = 0.25) {
  ss_end <- n_doses * tau_h
  et(amt  = dose_F,
     ii   = tau_h,
     addl = n_doses - 1,
     time = 0) %>%
    et(seq(0, ss_end + wash_h, by = dt))
}

# ── 4. Run Simulations & Extract Metrics ─────────────────────────────────────
# v2 핵심 변경: parent+metabolite 모델에서 Css_sim_parent_ngmL /
# Css_sim_met_ngmL / Css_sim_total_ngmL 세 개를 항상 별도로 저장.
# 더 이상 "Css_sim_ng" 하나에 Ctotal을 담아 parent인 것처럼 보이게 하지 않음.

results <- list()

# ── 4.1 Sertraline (1-cmt, 분리 불필요) ───────────────────────────────────────
p  <- pk$sertraline
ev <- make_events(p$F * p$dose_mg, p$tau_h, p$n_doses, p$wash_h)

sim_sert <- rxSolve(
  ode_1cmt, params = c(Ka = p$Ka, Vd = p$Vd, CL = p$CL),
  events = ev, inits = c(depot = 0, Cp = 0)
) %>% as.data.frame()

sim_sert_sd   <- sim_sert %>% filter(time < p$tau_h)
t_last_dose   <- (p$n_doses - 1) * p$tau_h
sim_sert_wash <- sim_sert %>% filter(time >= t_last_dose) %>%
  mutate(time = time - t_last_dose)

results[["Sertraline"]] <- list(
  Tmax_sim          = get_tmax_cmax(sim_sert_sd)$Tmax_h,
  thalf_sim         = get_thalf_sim(sim_sert_wash),
  Css_sim_parent_ng = round(get_css(
    sim_sert %>% filter(time >= t_last_dose, time <= p$n_doses * p$tau_h)) * 1000, 1),
  Css_sim_met_ng    = NA_real_,
  Css_sim_total_ng  = NA_real_,
  thalf_analytical  = round(log(2) * p$Vd / p$CL, 1)
)

# ── 4.2 Fluoxetine (parent+metabolite — 3-3 수정 적용) ───────────────────────
p  <- pk$fluoxetine
ev <- make_events(p$F * p$dose_mg, p$tau_h, p$n_doses, p$wash_h, dt = 1)

sim_flx <- rxSolve(
  ode_parent_met,
  params = c(Ka = p$Ka, Vd_p = p$Vd_p, CL_p = p$CL_p,
             Fm  = p$Fm, Vd_m = p$Vd_m, CL_m = p$CL_m),
  events = ev, inits = c(depot = 0, Cp = 0, Cm = 0)
) %>% as.data.frame()

sim_flx_sd   <- sim_flx %>% filter(time < p$tau_h)
t_last_dose  <- (p$n_doses - 1) * p$tau_h
sim_flx_wash <- sim_flx %>% filter(time >= t_last_dose) %>%
  mutate(time = time - t_last_dose)
sim_flx_ss   <- sim_flx %>% filter(time >= t_last_dose, time <= p$n_doses * p$tau_h)

css_flx_parent <- get_css(sim_flx_ss, "Cp")
css_flx_met    <- get_css(sim_flx_ss, "Cm")
css_flx_total  <- get_css(sim_flx_ss, "Ctotal")

results[["Fluoxetine"]] <- list(
  Tmax_sim          = get_tmax_cmax(sim_flx_sd, "Cp")$Tmax_h,
  thalf_sim         = get_thalf_sim(sim_flx_wash, "Cp"),
  Css_sim_parent_ng = round(css_flx_parent * 1000, 1),   # v2: parent만 (v1 버그 수정)
  Css_sim_met_ng    = round(css_flx_met    * 1000, 1),
  Css_sim_total_ng  = round(css_flx_total  * 1000, 1),   # 참고용으로 합계도 남김 (명시적으로 "total"이라 표기)
  thalf_analytical  = round(log(2) * p$Vd_p / p$CL_p, 1)
)

results[["Norfluoxetine"]] <- list(
  Tmax_sim          = NA,
  thalf_sim         = get_thalf_sim(sim_flx_wash, "Cm"),
  Css_sim_parent_ng = NA_real_,   # metabolite 자체 행이므로 parent/total 개념 없음
  Css_sim_met_ng    = round(css_flx_met * 1000, 1),
  Css_sim_total_ng  = NA_real_,
  thalf_analytical  = round(log(2) * p$Vd_m / p$CL_m, 1)
)

# ── 4.3 Paroxetine (1-cmt, 분리 불필요) ───────────────────────────────────────
p  <- pk$paroxetine
ev <- make_events(p$F * p$dose_mg, p$tau_h, p$n_doses, p$wash_h)

sim_par <- rxSolve(
  ode_1cmt, params = c(Ka = p$Ka, Vd = p$Vd, CL = p$CL),
  events = ev, inits = c(depot = 0, Cp = 0)
) %>% as.data.frame()

sim_par_sd   <- sim_par %>% filter(time < p$tau_h)
t_last_dose  <- (p$n_doses - 1) * p$tau_h
sim_par_wash <- sim_par %>% filter(time >= t_last_dose) %>%
  mutate(time = time - t_last_dose)

results[["Paroxetine"]] <- list(
  Tmax_sim          = get_tmax_cmax(sim_par_sd)$Tmax_h,
  thalf_sim         = get_thalf_sim(sim_par_wash),
  Css_sim_parent_ng = round(get_css(
    sim_par %>% filter(time >= t_last_dose, time <= p$n_doses * p$tau_h)) * 1000, 1),
  Css_sim_met_ng    = NA_real_,
  Css_sim_total_ng  = NA_real_,
  thalf_analytical  = round(log(2) * p$Vd / p$CL, 1)
)

# ── 4.4 Venlafaxine IR (parent+metabolite — 3-3 수정 적용) ───────────────────
p  <- pk$venlafaxine_IR
ev <- make_events(p$F * p$dose_mg, p$tau_h, p$n_doses, p$wash_h)

sim_ven_IR <- rxSolve(
  ode_parent_met,
  params = c(Ka = p$Ka, Vd_p = p$Vd_p, CL_p = p$CL_p,
             Fm  = p$Fm, Vd_m = p$Vd_m, CL_m = p$CL_m),
  events = ev, inits = c(depot = 0, Cp = 0, Cm = 0)
) %>% as.data.frame()

sim_ven_IR_sd   <- sim_ven_IR %>% filter(time < p$tau_h)
t_last_dose     <- (p$n_doses - 1) * p$tau_h
sim_ven_IR_wash <- sim_ven_IR %>% filter(time >= t_last_dose) %>%
  mutate(time = time - t_last_dose)
sim_ven_IR_ss   <- sim_ven_IR %>% filter(time >= t_last_dose, time <= p$n_doses * p$tau_h)

results[["Venlafaxine IR"]] <- list(
  Tmax_sim          = get_tmax_cmax(sim_ven_IR_sd, "Cp")$Tmax_h,
  thalf_sim         = get_thalf_sim(sim_ven_IR_wash, "Cp"),
  Css_sim_parent_ng = round(get_css(sim_ven_IR_ss, "Cp") * 1000, 1),  # v2: VEN parent만
  Css_sim_met_ng    = round(get_css(sim_ven_IR_ss, "Cm") * 1000, 1),  # ODV
  Css_sim_total_ng  = round(get_css(sim_ven_IR_ss, "Ctotal") * 1000, 1),
  thalf_analytical  = round(log(2) * p$Vd_p / p$CL_p, 1)
)

# ── 4.5 Venlafaxine XR (parent+metabolite — 3-3 수정 적용) ───────────────────
p  <- pk$venlafaxine_XR
ev <- make_events(p$F * p$dose_mg, p$tau_h, p$n_doses, p$wash_h)

sim_ven_XR <- rxSolve(
  ode_parent_met,
  params = c(Ka = p$Ka, Vd_p = p$Vd_p, CL_p = p$CL_p,
             Fm  = p$Fm, Vd_m = p$Vd_m, CL_m = p$CL_m),
  events = ev, inits = c(depot = 0, Cp = 0, Cm = 0)
) %>% as.data.frame()

sim_ven_XR_sd   <- sim_ven_XR %>% filter(time < p$tau_h)
t_last_dose     <- (p$n_doses - 1) * p$tau_h
sim_ven_XR_wash <- sim_ven_XR %>% filter(time >= t_last_dose) %>%
  mutate(time = time - t_last_dose)
sim_ven_XR_ss   <- sim_ven_XR %>% filter(time >= t_last_dose, time <= p$n_doses * p$tau_h)

results[["Venlafaxine XR"]] <- list(
  Tmax_sim          = get_tmax_cmax(sim_ven_XR_sd, "Cp")$Tmax_h,
  thalf_sim         = get_thalf_sim(sim_ven_XR_wash, "Cp"),
  Css_sim_parent_ng = round(get_css(sim_ven_XR_ss, "Cp") * 1000, 1),
  Css_sim_met_ng    = round(get_css(sim_ven_XR_ss, "Cm") * 1000, 1),
  Css_sim_total_ng  = round(get_css(sim_ven_XR_ss, "Ctotal") * 1000, 1),
  thalf_analytical  = round(log(2) * p$Vd_p / p$CL_p, 1)
)

results[["ODV"]] <- list(
  Tmax_sim          = NA,
  thalf_sim         = get_thalf_sim(sim_ven_XR_wash, "Cm"),
  Css_sim_parent_ng = NA_real_,
  Css_sim_met_ng    = round(get_css(sim_ven_XR_ss, "Cm") * 1000, 1),
  Css_sim_total_ng  = NA_real_,
  thalf_analytical  = round(log(2) * pk$venlafaxine_XR$Vd_m /
                             pk$venlafaxine_XR$CL_m, 1)
)

# ── 5. parent:metabolite 비율 (3-3 수정의 핵심 산출물) ────────────────────────
# 원고 Table 2 / Section 3.1 / Limitation 7 / 3.9절에 쓰인 "NFX:FLX = 0.51"이
# 실제로는 (met)/(parent+met)에 가까운 잘못된 계산이었음을 바로잡기 위한 표.
# true ratio = met / parent (예: NFX/FLX, ODV/VEN).

ratio_table <- tibble(
  pair = c("Norfluoxetine : Fluoxetine (parent)",
           "ODV : Venlafaxine IR (parent)",
           "ODV : Venlafaxine XR (parent)"),
  met_ngmL = c(
    results[["Fluoxetine"]]$Css_sim_met_ng,
    results[["Venlafaxine IR"]]$Css_sim_met_ng,
    results[["Venlafaxine XR"]]$Css_sim_met_ng
  ),
  parent_ngmL = c(
    results[["Fluoxetine"]]$Css_sim_parent_ng,
    results[["Venlafaxine IR"]]$Css_sim_parent_ng,
    results[["Venlafaxine XR"]]$Css_sim_parent_ng
  ),
  total_ngmL = c(
    results[["Fluoxetine"]]$Css_sim_total_ng,
    results[["Venlafaxine IR"]]$Css_sim_total_ng,
    results[["Venlafaxine XR"]]$Css_sim_total_ng
  )
) %>%
  mutate(
    true_ratio_met_over_parent  = round(met_ngmL / parent_ngmL, 3),
    wrong_ratio_v1_met_over_sum = round(met_ngmL / total_ngmL, 3)  # v1이 실수로 보고했던 방식 (대조용)
  )

ratio_path <- file.path(OUT_DIR, "parent_metabolite_ratio_v2.csv")
write_csv(ratio_table, ratio_path)
cat("\n── Parent:Metabolite ratio check (3-3 fix) ──\n")
print(ratio_table)
cat(sprintf("✓ Saved: %s\n", ratio_path))

# ── 6. Label Reference Table (parent 기준으로 통일) ──────────────────────────
# v1에서 Fluoxetine/Venlafaxine의 Css_lo/Css_hi가 "Ctotal 추정" 범위였던 것을
# v2에서는 가능한 한 parent 단독 라벨 범위로 맞춤 (불확실한 항목은 "—"로 비움
# — 저자가 실제 라벨 수치로 재확인해야 함, SOURCE_CHECKS.md 참고).

label <- tribble(
  ~drug,             ~std_dose_str, ~Tmax_lo, ~Tmax_hi,
  ~thalf_lo, ~thalf_hi,
  ~Css_lo,   ~Css_hi,    ~label_source,

  "Sertraline",      "50 mg QD",    4.5,  8.4,   22,   36,    20,   60,
  "FDA label (Zoloft); Vd 출처 재확인 필요(unsourced_pending)",

  "Fluoxetine",      "20 mg QD",    6.0,  8.0,  120,  192,    NA,   NA,
  "FDA label (Prozac); parent 단독 Css 범위는 저자 재확인 필요 — v1의 Ctotal 추정 범위 제거",

  "Norfluoxetine",   "—",           NA,   NA,   192,  384,    NA,   NA,
  "FDA label (Prozac)",

  "Paroxetine",      "20 mg QD",    3.0,  8.0,    7,   65,    30,  100,
  "FDA label (Paxil); Li 2022 (Tmax only)",

  "Venlafaxine IR",  "25 mg TID",   1.0,  2.5,    3,    7,    NA,   NA,
  "FDA label (Effexor XR 2012); parent 단독 Css 범위는 저자 재확인 필요",

  "Venlafaxine XR",  "75 mg QD",    4.0,  8.0,    3,    7,    NA,   NA,
  "FDA label (Effexor XR 2012); parent 단독 Css 범위는 저자 재확인 필요",

  "ODV",             "—",           NA,   NA,     9,   13,    NA,   NA,
  "FDA label (Effexor XR 2012)"
)

# ── 7. Assemble Qualification Table (v2) ─────────────────────────────────────

drug_order <- c("Sertraline", "Fluoxetine", "Norfluoxetine",
                "Paroxetine", "Venlafaxine IR", "Venlafaxine XR", "ODV")

qual_rows <- map_dfr(drug_order, function(d) {
  r   <- results[[d]]
  lab <- label %>% filter(drug == d)

  flag_thalf <- case_when(
    is.na(r$thalf_sim)                              ~ "N/A",
    r$thalf_sim >= lab$thalf_lo &
      r$thalf_sim <= lab$thalf_hi                   ~ "PASS",
    r$thalf_sim  < lab$thalf_lo * 0.75 |
      r$thalf_sim  > lab$thalf_hi * 1.25            ~ "CONCERN",
    TRUE                                            ~ "BORDERLINE"
  )

  flag_tmax <- case_when(
    is.na(r$Tmax_sim) | is.na(lab$Tmax_lo)         ~ "N/A",
    r$Tmax_sim >= lab$Tmax_lo &
      r$Tmax_sim <= lab$Tmax_hi                     ~ "PASS",
    TRUE                                            ~ "BORDERLINE"
  )

  # v2: parent 기준 Css로 flag 계산 (분리된 값 사용, total 아님)
  css_ref <- if (!is.na(r$Css_sim_parent_ng)) r$Css_sim_parent_ng else r$Css_sim_met_ng

  flag_css <- case_when(
    is.na(css_ref) | is.na(lab$Css_lo)             ~ "N/A",
    css_ref >= lab$Css_lo & css_ref <= lab$Css_hi   ~ "PASS",
    css_ref < lab$Css_lo * 0.5 |
      css_ref > lab$Css_hi * 2                      ~ "CONCERN",
    TRUE                                            ~ "BORDERLINE"
  )

  overall <- case_when(
    flag_thalf == "CONCERN" | flag_css == "CONCERN" ~ "CONCERN — review",
    flag_thalf == "PASS"                            ~ "PASS",
    flag_thalf == "BORDERLINE"                      ~ "BORDERLINE",
    TRUE                                            ~ "N/A"
  )

  tibble(
    Drug                  = d,
    Std_dose              = lab$std_dose_str,
    Tmax_sim_h            = r$Tmax_sim,
    Tmax_label_range_h    = ifelse(is.na(lab$Tmax_lo), "—",
                                   paste0(lab$Tmax_lo, "–", lab$Tmax_hi)),
    Tmax_flag             = flag_tmax,
    thalf_sim_h           = r$thalf_sim,
    thalf_analytical_h    = r$thalf_analytical,
    thalf_label_range_h   = paste0(lab$thalf_lo, "–", lab$thalf_hi),
    thalf_flag            = flag_thalf,
    Css_sim_parent_ngmL   = r$Css_sim_parent_ng,   # v2: 항상 parent 단독
    Css_sim_met_ngmL      = r$Css_sim_met_ng,      # v2: 항상 metabolite 단독
    Css_sim_total_ngmL    = r$Css_sim_total_ng,    # v2: 참고용 합계 (명시적 라벨)
    Css_label_range_ngmL  = ifelse(is.na(lab$Css_lo), "—",
                                   paste0(lab$Css_lo, "–", lab$Css_hi)),
    Css_flag              = flag_css,
    Overall               = overall,
    Source                = lab$label_source
  )
})

# ── 8. Save CSV ───────────────────────────────────────────────────────────────
csv_path <- file.path(OUT_DIR, "model_qualification_table_v2.csv")
write_csv(qual_rows, csv_path)
message("✓ Qualification table (v2) saved: ", csv_path)

# ── 9. Print to console ───────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════\n")
cat(" MODEL QUALIFICATION TABLE v2  (Phase 3.5, 3-3 fix applied)\n")
cat("═══════════════════════════════════════════════════════\n\n")

print_cols <- c("Drug", "thalf_sim_h", "thalf_analytical_h",
                "thalf_label_range_h", "thalf_flag",
                "Css_sim_parent_ngmL", "Css_sim_met_ngmL", "Css_sim_total_ngmL",
                "Overall")
print(qual_rows[, print_cols], n = Inf, width = Inf)

# ── 10. Visualisation (동일 로직, v1과 구조 동일) ────────────────────────────

plot_data <- qual_rows %>%
  filter(!is.na(thalf_sim_h)) %>%
  left_join(label %>% select(drug, thalf_lo, thalf_hi),
            by = c("Drug" = "drug")) %>%
  mutate(Drug = factor(Drug, levels = rev(drug_order)))

p_thalf <- ggplot(plot_data, aes(y = Drug)) +
  geom_rect(aes(xmin = thalf_lo, xmax = thalf_hi,
                ymin = as.numeric(Drug) - 0.4,
                ymax = as.numeric(Drug) + 0.4),
            fill = "#B8D4E8", alpha = 0.6) +
  geom_point(aes(x = thalf_sim_h, colour = thalf_flag),
             size = 4, shape = 18) +
  geom_point(aes(x = thalf_analytical_h),
             size = 3, shape = 5, colour = "grey40") +
  scale_colour_manual(
    values = c("PASS" = "#2E8B57", "BORDERLINE" = "#DAA520",
               "CONCERN" = "#C0392B", "N/A" = "grey60"),
    name = "t½ flag"
  ) +
  scale_x_log10(
    breaks = c(1, 3, 5, 10, 20, 50, 100, 200, 400),
    labels = scales::comma
  ) +
  labs(
    title    = "Model Qualification v2: Simulated t½ vs FDA Label Range",
    subtitle = "Blue band = label range  |  ◆ = simulated  |  ⬦ = analytical (0.693·Vd/CL)  |  F=1 (D3), parent/met separated (3-3 fix)",
    x        = "Terminal half-life (h, log scale)",
    y        = NULL,
    caption  = paste0("Phase 3.5 qualification v2  |  seed = 2025  |  ",
                      format(Sys.Date(), "%Y-%m-%d"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

plot_path <- file.path(OUT_DIR, "model_qualification_plot_v2.png")
ggsave(plot_path, p_thalf, width = 10, height = 6, dpi = 300)
message("✓ Qualification plot (v2) saved: ", plot_path)

# ── 11. Flag summary ──────────────────────────────────────────────────────────
cat("\n── Overall flags ──\n")
qual_rows %>%
  select(Drug, thalf_flag, Tmax_flag, Css_flag, Overall) %>%
  print(n = Inf)

n_concern <- sum(qual_rows$Overall == "CONCERN — review", na.rm = TRUE)
if (n_concern == 0) {
  cat("\n✓ All drugs PASS or BORDERLINE. Safe to proceed to v2 Monte Carlo.\n\n")
} else {
  cat(sprintf(
    "\n⚠ %d drug(s) flagged CONCERN. Review parameters before proceeding.\n\n",
    n_concern
  ))
}

cat("▶ 다음 단계: mc_run_v2.R 실행 (N=100 테스트 → N=10000 본실행), 그다음 sensitivity_v2.R\n")
