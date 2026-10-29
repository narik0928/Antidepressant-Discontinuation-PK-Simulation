# =============================================================================
# mc_run_v3.R -- Phase 4 Monte Carlo Runner, v3
# Author: Nari Kim (code assist: Claude)
#
# HANDOFF_TO_NARI_v3.md 7-1-2 요구사항 반영. mc_params_v3.R(z값 역산 +
# D11 중앙값보존 재샘플링 + CYP2D6 재배정)의 출력을 읽어 ODE를 풀고 지표를
# 뽑는다 -- mc_run_v2.R이 mc_params_all_n10000_v1.csv를 읽던 것과 같은 구조.
#
# v2 대비 바뀐 것
#   D6  : paroxetine 임상 투여를 20mg -> 30mg/day로 올림 (censoring 감소
#         목적, HANDOFF_v3 6절).
#   B6  : paroxetine 관찰창을 14일 -> 30일로 늘림 (censoring 감소 목적).
#   6절 : sertraline 관찰창을 30일 -> 14일로 좁힘 (HANDOFF_v3 6절 표 --
#         시작점이며, t25% censoring이 1%를 넘으면 늘리고 보고할 것).
#   D5  : CYP2D6 보정은 mc_params_v3.R 단계에서 이미 CL_sampled에
#         반영돼 있음(paroxetine PM x0.58만, 나머지는 배수 없음) --
#         이 스크립트는 phenotype 배수를 다시 적용하지 않는다.
#   D11 : 파라미터 자체가 이미 중앙값보존 샘플링으로 재적용돼 있음(이
#         스크립트는 소비만 함, 변환하지 않음).
#   정밀도: C0, AUC(절대/정규화), frac96/frac120은 반올림하지 않고 원래
#         정밀도 그대로 저장(HANDOFF_v3 7-1-2: "유효숫자 6자리 이상").
#         v2는 round(..., 6)을 썼지만 정규화값 역산이 막혔던 문제가
#         있었음(7-1-2) -- v3는 애초에 반올림을 하지 않는다.
#   신규 지표: 96h/120h 잔존 비율 C(t)/C0 (frac96, frac120) 추가.
#         active-metabolite 약물(fluoxetine, venlafaxine)은 total active
#         moiety(Ctotal = Cp + Cm) 기준 -- v2와 동일.
#   삭제: legacy_Css, reldiff_72h/pct_drop_72h -- v1-v2 비교용(S1)이었고
#         v3 7-1-2가 요구하는 지표 목록에 없어 제거함. 필요하면 원시 sim
#         출력에서 다시 계산 가능.
#
# Input : outputs/mc_v3/mc_params_v3_all_n<N>.csv (mc_params_v3.R 산출물)
# Output: outputs/mc_v3/mc_results_v3_<drug>_n<N>.csv
#         outputs/mc_v3/mc_results_v3_all_n<N>.csv
#
# 실행: source("R/simulation/mc_run_v3.R")
# 먼저 N=100(mc_params_v3.R의 N을 100으로 바꿔 재실행 후)으로 에러 없이
# 끝까지 도는지, 그리고 calibration/typical_v3.py 기대값과 5% 이내로
# 맞는지 확인한 뒤 N=10000 본 실행. (HANDOFF_v3 7-2 -- 5% 넘게 다르면
# 멈추고 덕현과 공유할 것)
# =============================================================================

suppressPackageStartupMessages({
  library(rxode2)
  library(tidyverse)
})

# ── 0. Setup ──────────────────────────────────────────────────────────────────
BASE     <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR  <- file.path(BASE, "outputs", "mc_v3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N <- 100    # mc_params_v3.R에서 생성한 N과 일치해야 함. 먼저 100으로 점검.

source(file.path(BASE, "data", "pk_parameters", "pk_parameters_v3.R"))

# ── 1. ODE Models (v2와 동일 구조, 변경 없음) ─────────────────────────────────

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

# ── 2. Drug Configuration (v3, HANDOFF_v3 6절) ────────────────────────────────
# F는 pk_params_v3에서 읽음(전부 1.0, D3 유지). CYP2D6 배수는 mc_params_v3.R
# 단계에서 이미 CL_sampled에 반영돼 있으므로 여기서는 다시 곱하지 않는다.

drug_config <- list(

  Sertraline = list(
    model    = "1cmt",
    Ka       = pk_params$sertraline$Ka,
    F        = pk_params$sertraline$F,
    Fm       = NA,
    dose_mg  = 50,
    tau_h    = 24,
    n_doses  = 14,
    obs_h    = 30 * 24,   # 추가 지시 1(2-2) -- v3 본문 6절의 "14일"은 덕현 쪽
                           # 착오로 확인됨. v2와 같은 30일로 되돌림(14일이면
                           # t10%에서 약 50명 censoring, 30일이면 0명).
    dt       = 1.0,
    conc_col = "Cp"
  ),

  Fluoxetine = list(
    model    = "parent_met",
    Ka       = pk_params$fluoxetine$Ka,
    F        = pk_params$fluoxetine$F,
    Fm       = pk_params$fluoxetine$norfluoxetine$Fm,
    dose_mg  = 20,
    tau_h    = 24,
    n_doses  = 60,
    obs_h    = 400 * 24,  # v3 6절: 관찰창 약 400일 (D4, 유지)
    dt       = 2.0,
    conc_col = "Ctotal"
  ),

  Paroxetine = list(
    model    = "1cmt",
    Ka       = pk_params$paroxetine$Ka,
    F        = pk_params$paroxetine$F,
    Fm       = NA,
    dose_mg  = 30,        # v3 D6: 20mg -> 30mg/day (라벨 정상상태 PK 자료 용량)
    tau_h    = 24,
    n_doses  = 14,
    obs_h    = 30 * 24,   # v3 6절/B6: 관찰창 14일 -> 30일 (censoring 감소 목적)
    dt       = 0.5,
    conc_col = "Cp"
  ),

  Venlafaxine_IR = list(
    model    = "parent_met",
    Ka       = pk_params$venlafaxine_IR$Ka,
    F        = pk_params$venlafaxine_IR$F,
    Fm       = pk_params$venlafaxine_IR$ODV$Fm,
    dose_mg  = 25,     # 25mg TID = 75mg/day (독립 점검용, D7)
    tau_h    = 8,
    n_doses  = 21,
    obs_h    = 14 * 24,   # 추가 지시 1(2-2) -- 7일 -> 14일. v3 파라미터로는
                           # 반감기가 길어져(parent 5.5h, ODV 13.1h) 7일이면
                           # t10%에서 IR 32명 censoring(t25%에서도 2명),
                           # 14일이면 1명. AUC deficit 적분구간(168h)은 그대로.
    dt       = 0.25,
    conc_col = "Ctotal"
  ),

  Venlafaxine_XR = list(
    model    = "parent_met",
    Ka       = pk_params$venlafaxine_XR$Ka,
    F        = pk_params$venlafaxine_XR$F,
    Fm       = pk_params$venlafaxine_XR$ODV$Fm,
    dose_mg  = 75,
    tau_h    = 24,
    n_doses  = 7,
    obs_h    = 14 * 24,   # 추가 지시 1(2-2) -- 7일 -> 14일. 7일이면 t10%에서
                           # XR 23명 censoring, 14일이면 0명.
    dt       = 0.25,
    conc_col = "Ctotal"
  )
)

# ── 3. Helper functions ───────────────────────────────────────────────────────

find_threshold_time_v3 <- function(times, concs, threshold) {
  if (length(times) < 2 || is.na(threshold) || threshold <= 0) {
    return(list(t = NA_real_, status = "calc_failed"))
  }
  idx <- which(concs <= threshold)[1]
  if (is.na(idx)) {
    return(list(t = NA_real_, status = "censored"))        # 관찰 종료까지 미도달 (우측 검열)
  }
  if (idx == 1) {
    return(list(t = NA_real_, status = "calc_failed"))      # 이론상 발생하면 안 되는 케이스
  }
  t1 <- times[idx - 1]; t2 <- times[idx]
  c1 <- concs[idx - 1]; c2 <- concs[idx]
  list(t = t1 + (threshold - c1) * (t2 - t1) / (c2 - c1), status = "reached")
}

trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# 환자 1명 PK metrics 추출 (v3)
# C0: 중단 시점 순간 농도 (D1, 유지)
extract_patient_metrics_v3 <- function(.x, C0, t_half_impl_parent, t_half_impl_met) {

  times <- .x$post_time
  concs <- .x$conc

  base_na <- tibble(
    C0 = NA_real_,
    reldiff_24h = NA_real_,
    frac96 = NA_real_, frac120 = NA_real_,
    t_75pct = NA_real_, status_75pct = "calc_failed",
    t_50pct = NA_real_, status_50pct = "calc_failed",
    t_25pct = NA_real_, status_25pct = "calc_failed",
    t_10pct = NA_real_, status_10pct = "calc_failed",
    AUC_deficit_7d = NA_real_, AUC_deficit_7d_norm_pct = NA_real_,
    t_half_impl_parent = t_half_impl_parent, t_half_impl_met = t_half_impl_met
  )

  if (is.na(C0) || C0 <= 0 || length(times) < 3) return(base_na)

  C_24h  <- approx(times, concs, xout = 24)$y
  C_96h  <- approx(times, concs, xout = 96)$y
  C_120h <- approx(times, concs, xout = 120)$y

  # 24h 평균 상대변화율 (HANDOFF_v3 7-1-2: "24 h 평균 상대변화율", v2와
  # 동일 정의 -- (C(24h)-C0)/(24*C0))
  reldiff_24h <- if (!is.na(C_24h)) (C_24h - C0) / (24 * C0) else NA_real_

  # 96h/120h 잔존 비율 C(t)/C0 (HANDOFF_v3 7-2 기대값 표와 비교하는 지표)
  frac96  <- if (!is.na(C_96h))  C_96h  / C0 else NA_real_
  frac120 <- if (!is.na(C_120h)) C_120h / C0 else NA_real_

  r75 <- find_threshold_time_v3(times, concs, 0.75 * C0)
  r50 <- find_threshold_time_v3(times, concs, 0.50 * C0)
  r25 <- find_threshold_time_v3(times, concs, 0.25 * C0)
  r10 <- find_threshold_time_v3(times, concs, 0.10 * C0)

  # AUC deficit -- t50% 도달 여부와 무관하게 항상 계산 (7일 또는 obs_h 중 짧은 쪽)
  obs_max <- min(168, max(times, na.rm = TRUE))
  post_7d <- .x %>% filter(post_time <= obs_max)
  if (nrow(post_7d) >= 2) {
    AUC_actual     <- trapz(post_7d$post_time, post_7d$conc)
    AUC_expected   <- C0 * obs_max
    AUC_deficit_7d <- AUC_expected - AUC_actual
    AUC_deficit_7d_norm_pct <- 100 * AUC_deficit_7d / AUC_expected
  } else {
    AUC_deficit_7d <- NA_real_
    AUC_deficit_7d_norm_pct <- NA_real_
  }

  # 정밀도: C0/AUC/frac은 반올림하지 않음 (HANDOFF_v3 7-1-2, 유효숫자 6자리 이상).
  # 시간(h)·반감기만 6자리로 반올림 (round-trip에 지장 없는 보조 정보).
  tibble(
    C0                       = C0,
    reldiff_24h              = reldiff_24h,
    frac96                   = frac96,
    frac120                  = frac120,
    t_75pct                  = round(r75$t, 6), status_75pct = r75$status,
    t_50pct                  = round(r50$t, 6), status_50pct = r50$status,
    t_25pct                  = round(r25$t, 6), status_25pct = r25$status,
    t_10pct                  = round(r10$t, 6), status_10pct = r10$status,
    AUC_deficit_7d           = AUC_deficit_7d,
    AUC_deficit_7d_norm_pct  = AUC_deficit_7d_norm_pct,
    t_half_impl_parent       = round(t_half_impl_parent, 6),
    t_half_impl_met          = if (!is.na(t_half_impl_met)) round(t_half_impl_met, 6) else NA_real_
  )
}

# ── 4. 배치 시뮬레이션 함수 (v3) ──────────────────────────────────────────────

run_batch_v3 <- function(drug_name, batch_df, cfg) {

  t_disc <- cfg$n_doses * cfg$tau_h

  ev <- et(amt  = cfg$F * cfg$dose_mg,
           ii   = cfg$tau_h,
           addl = cfg$n_doses - 1,
           time = 0) %>%
    et(seq(0, t_disc + cfg$obs_h, by = cfg$dt))

  if (cfg$model == "1cmt") {
    model <- ode_1cmt
    inits <- c(depot = 0, Cp = 0)
    get_params <- function(row) c(Ka = cfg$Ka, Vd = row$Vd_sampled, CL = row$CL_sampled)
  } else {
    model <- ode_parent_met
    inits <- c(depot = 0, Cp = 0, Cm = 0)
    get_params <- function(row) c(Ka = cfg$Ka, Vd_p = row$Vd_sampled, CL_p = row$CL_sampled,
                                   Fm = cfg$Fm, Vd_m = row$Vd_met_sampled, CL_m = row$CL_met_sampled)
  }

  sim_list <- lapply(seq_len(nrow(batch_df)), function(i) {
    row   <- batch_df[i, , drop = FALSE]
    sim_i <- rxSolve(model, get_params(row), ev, inits = inits) %>% as.data.frame()
    sim_i$id <- row$id
    sim_i
  })

  sim <- bind_rows(sim_list) %>%
    rename_with(~ "conc", .cols = all_of(cfg$conc_col))

  # C0: 중단 시점(t_disc)과 가장 가까운 격자점의 농도 (D1, 유지)
  c0_df <- sim %>%
    group_by(id) %>%
    slice(which.min(abs(time - t_disc))) %>%
    ungroup() %>%
    transmute(id, C0 = conc)

  post_df <- sim %>%
    filter(time >= t_disc) %>%
    mutate(post_time = time - t_disc) %>%
    select(id, post_time, conc) %>%
    left_join(c0_df, by = "id")

  # 환자별 구현 반감기 (D11 확인용 -- 7-2 마지막 문단)
  half_life_df <- batch_df %>%
    transmute(
      id,
      t_half_impl_parent = 0.693 * Vd_sampled / CL_sampled,
      t_half_impl_met = if (cfg$model == "parent_met") 0.693 * Vd_met_sampled / CL_met_sampled else NA_real_
    )

  metrics <- post_df %>%
    left_join(half_life_df, by = "id") %>%
    group_by(id) %>%
    group_modify(~ extract_patient_metrics_v3(
      .x,
      C0 = unique(.x$C0),
      t_half_impl_parent = unique(.x$t_half_impl_parent),
      t_half_impl_met = unique(.x$t_half_impl_met)
    )) %>%
    ungroup()

  batch_df %>%
    select(id, drug, cyp2d6, cyp2d6_v1_label, CL_sampled, Vd_sampled, sampling_mode) %>%
    left_join(metrics, by = "id")
}

# ── 5. 전체 실행 ──────────────────────────────────────────────────────────────
# mc_params_v3.R의 출력을 읽는다 (새로 샘플링하지 않음 -- mc_params_v3.R이
# z값 역산/재적용을 이미 끝낸 상태).

params_path <- file.path(BASE, "outputs", "mc_v3", sprintf("mc_params_v3_all_n%d.csv", N))
if (!file.exists(params_path)) {
  stop(sprintf(
    "%s 없음. 먼저 R/simulation/mc_params_v3.R을 N=%d로 실행해 이 파일을 만드세요.",
    params_path, N
  ))
}
params_all <- read_csv(params_path, show_col_types = FALSE)

BATCH_SIZE <- 1000
results_list <- list()

for (drug_name in names(drug_config)) {

  cfg         <- drug_config[[drug_name]]
  drug_params <- params_all %>% filter(drug == drug_name)
  n_batches   <- ceiling(nrow(drug_params) / BATCH_SIZE)

  cat(sprintf("\n[v3] %s: %d patients / %d batches (dose=%gmg, obs_h=%g h, dt=%g h)\n",
              drug_name, nrow(drug_params), n_batches, cfg$dose_mg, cfg$obs_h, cfg$dt))

  t_start      <- proc.time()["elapsed"]
  drug_results <- vector("list", n_batches)

  for (b in seq_len(n_batches)) {
    idx_start <- (b - 1) * BATCH_SIZE + 1
    idx_end   <- min(b * BATCH_SIZE, nrow(drug_params))
    batch_df  <- drug_params[idx_start:idx_end, ]
    drug_results[[b]] <- run_batch_v3(drug_name, batch_df, cfg)
    cat(sprintf("  batch %02d/%02d done (n=%d)\n", b, n_batches, nrow(batch_df)))
  }

  drug_df <- bind_rows(drug_results)
  elapsed <- round(proc.time()["elapsed"] - t_start, 1)

  out_path <- file.path(OUT_DIR, sprintf("mc_results_v3_%s_n%d.csv",
                                          tolower(gsub("_", "", drug_name)), N))
  write_csv(drug_df, out_path)
  cat(sprintf("Saved: %s  [%.1f sec]\n", out_path, elapsed))

  results_list[[drug_name]] <- drug_df
}

# ── 6. 전체 합본 저장 ─────────────────────────────────────────────────────────

results_all <- bind_rows(results_list)
all_path    <- file.path(OUT_DIR, sprintf("mc_results_v3_all_n%d.csv", N))
write_csv(results_all, all_path)

cat(sprintf("\nmc_run_v3.R 완료\n"))
cat(sprintf("   Total rows : %d\n", nrow(results_all)))
cat(sprintf("   Saved      : %s\n", all_path))

# ── 7. 빠른 검증 -- status 코드별 분포 + 지표 요약 + 7-2 기대값 대조 ─────────

cat("\n-- Status code distribution (t_50pct) --\n")
results_all %>%
  count(drug, status_50pct) %>%
  group_by(drug) %>%
  mutate(pct = round(100 * n / sum(n), 2)) %>%
  ungroup() %>%
  print(n = Inf)

cat("\n-- t25% censoring rate (HANDOFF_v3 6절: 1% 넘으면 관찰창을 늘리고 보고) --\n")
results_all %>%
  group_by(drug) %>%
  summarise(pct_censored_25 = round(100 * mean(status_25pct == "censored", na.rm = TRUE), 2), .groups = "drop") %>%
  print(n = Inf)

cat("\n-- Quick validation (median metrics by drug, reached 환자만) --\n")
results_all %>%
  filter(status_50pct == "reached") %>%
  group_by(drug) %>%
  summarise(
    C0_ngmL         = round(median(C0, na.rm = TRUE) * 1000, 2),
    t50_h           = round(median(t_50pct, na.rm = TRUE), 2),
    n_reached       = n(),
    AUCdef_norm_pct = round(median(AUC_deficit_7d_norm_pct, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

cat("\n-- D11 재확인: 환자별 구현 반감기 중앙값 (sertraline 전체, paroxetine normal만) --\n")
cat("기대값: sertraline 27.2h, paroxetine(normal) 21.0h -- 2% 이내여야 함\n")
results_all %>%
  filter((drug == "Sertraline") | (drug == "Paroxetine" & cyp2d6 == "normal")) %>%
  group_by(drug) %>%
  summarise(t_half_impl_median = round(median(t_half_impl_parent, na.rm = TRUE), 3), .groups = "drop") %>%
  print()

cat("\n다음 단계: N=100으로 이 출력을 calibration/typical_v3.py 기대값(7-2)과 5% 이내로\n")
cat("비교 -- 맞으면 mc_params_v3.R의 N을 10000으로 바꿔 본 실행, 그다음\n")
cat("qualification_v3.R, sensitivity_v3.R, stability_v3.R, tables_v3.R\n\n")
