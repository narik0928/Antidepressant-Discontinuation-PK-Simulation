# =============================================================================
# mc_run_v2.R — Phase 4 Monte Carlo Runner, v2
# Author: Nari Kim (코드 작성 보조: Claude)
#
# 원본 R/simulation/mc_run_single.R은 건드리지 않음. 새 파일로 분리.
# 덕현 인계 문서(HANDOFF_TO_NARI.md 6-3, RESUBMISSION_PLAN.md 3단계)의
# D1·D3·D4 확정 사항을 반영한 재실행 스크립트.
#
# v1 대비 바뀐 것
#   D1 (a): C0를 "중단 시점의 실제(순간) 농도"로 계산. 기존 코드는 마지막
#           투약 간격의 평균 농도(Css)를 기준으로 썼음 — 원고 2.5절의
#           정의("plasma concentration at discontinuation")와 불일치했던
#           부분. 이번엔 시뮬레이션 격자가 t_disc와 정확히 맞아떨어지는
#           시점의 농도값을 그대로 C0로 씀. 기존 방식의 평균 Css도
#           legacy_Css 컬럼으로 같이 저장함 (S1: 구-신 정의 짝지은 비교용).
#   D3    : 모든 약물 F=1 (pk_parameters_v2.R에서 통일, 이 스크립트는
#           F 값을 따로 적지 않고 cfg$F로 읽음).
#   D4(i) : fluoxetine 파라미터는 그대로 두고, 관찰창을 60일→약 400일로
#           늘리고 dt를 1h→2h로 넓힘(계산량 때문). "60일=steady state"
#           주장은 하지 않고, 도달 못한 환자는 censored로 보고함.
#   공통  : <1h 제외 규칙 없음. 환자별 status 코드(reached/censored/
#           calc_failed) 저장. t50%/t25%/24h diff/절대 AUC deficit/
#           정규화 AUC deficit을 전부 같은 C0 기준, 환자별로 계산.
#           t50% 미도달이어도 AUC는 계산함. 환자별 구현 반감기도 저장.
#
# Input : outputs/mc/mc_params_all_n10000_v1.csv (6-1단계에서 백업한
#         기존 가상환자 코호트 — 새로 샘플링하지 않음. 정의 변경 효과와
#         표본 차이를 분리하기 위함, GPT_STAGE2_COMPARISON.md 5절 권고1)
# Output: outputs/mc_v2/mc_results_v2_<drug>_n<N>.csv
#         outputs/mc_v2/mc_results_v2_all_n<N>.csv
#
# 실행: source("R/simulation/mc_run_v2.R")
# 먼저 N을 100으로 낮춰서 에러 없이 끝까지 도는지 확인한 뒤 10000으로.
# =============================================================================

suppressPackageStartupMessages({
  library(rxode2)
  library(tidyverse)
})

# ── 0. Setup ──────────────────────────────────────────────────────────────────
BASE     <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR  <- file.path(BASE, "outputs", "mc_v2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(2025)

N          <- 100     # 먼저 100으로 테스트. 통과하면 10000으로 바꿔서 본 실행.
BATCH_SIZE <- 1000

source(file.path(BASE, "data", "pk_parameters", "pk_parameters_v2.R"))

# ── 1. ODE Models (v1과 동일 구조, 변경 없음) ─────────────────────────────────

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

# ── 2. Drug Configuration (v2) ────────────────────────────────────────────────
# F는 pk_params에서 그대로 읽음 (v2: 전부 1.0). dose_mg·tau_h·n_doses는
# D2에 따라 숫자 유지. Fluoxetine만 obs_h·dt가 D4(i)에 따라 바뀜.

drug_config <- list(

  Sertraline = list(
    model    = "1cmt",
    Ka       = pk_params$sertraline$Ka,
    F        = pk_params$sertraline$F,          # v2: 1.0 (변경 없음)
    Fm       = NA,
    dose_mg  = 50,
    tau_h    = 24,
    n_doses  = 14,
    obs_h    = 30 * 24,
    dt       = 1.0,
    conc_col = "Cp"
  ),

  Fluoxetine = list(
    model    = "parent_met",
    Ka       = pk_params$fluoxetine$Ka,
    F        = pk_params$fluoxetine$F,                     # v2: 0.70 → 1.0
    Fm       = pk_params$fluoxetine$norfluoxetine$Fm,
    dose_mg  = 20,
    tau_h    = 24,
    n_doses  = 60,
    obs_h    = 400 * 24,   # v2 (D4-i): 60일 → 약 400일. PM 99백분위
                            # t25%(RESUBMISSION_PLAN.md 0-1절: 99백분위
                            # 반감기 192일, 최댓값 372일)를 덮도록 설정.
    dt       = 2.0,         # v2: 1h → 2h. 관찰창이 6.7배 늘어난 만큼
                            # 계산량 조절. t50%/t25% 보간 정확도에는
                            # 큰 영향 없을 것으로 예상(감소 시간축이 일 단위).
    conc_col = "Ctotal"
  ),

  Paroxetine = list(
    model    = "1cmt",
    Ka       = pk_params$paroxetine$Ka,
    F        = pk_params$paroxetine$F,          # v2: 0.50 → 1.0
    Fm       = NA,
    dose_mg  = 20,
    tau_h    = 24,
    n_doses  = 14,
    obs_h    = 14 * 24,
    dt       = 0.5,
    conc_col = "Cp"
  ),

  Venlafaxine_IR = list(
    model    = "parent_met",
    Ka       = pk_params$venlafaxine_IR$Ka,
    F        = pk_params$venlafaxine_IR$F,      # v2: 0.45 → 1.0
    Fm       = pk_params$venlafaxine_IR$ODV$Fm,
    dose_mg  = 25,     # 25mg TID = 75mg/day
    tau_h    = 8,
    n_doses  = 21,
    obs_h    = 7 * 24,
    dt       = 0.25,
    conc_col = "Ctotal"
  ),

  Venlafaxine_XR = list(
    model    = "parent_met",
    Ka       = pk_params$venlafaxine_XR$Ka,
    F        = pk_params$venlafaxine_XR$F,      # v2: 0.45 → 1.0
    Fm       = pk_params$venlafaxine_XR$ODV$Fm,
    dose_mg  = 75,
    tau_h    = 24,
    n_doses  = 7,
    obs_h    = 7 * 24,
    dt       = 0.25,
    conc_col = "Ctotal"
  )
)

# ── 3. Helper functions ───────────────────────────────────────────────────────

# 선형 보간으로 threshold 도달 시간 계산 + status 코드.
# v1과 차이: idx==1(관찰 시작부터 이미 기준 이하)을 더 이상 "즉시 도달"로
# 조용히 처리하지 않음 — v2에서는 post_time=0의 농도가 정의상 C0 자신이라
# 이론적으로 idx==1이 나오면 안 됨(1보다 작을 수 없으므로). 혹시 수치
# 이슈로 발생하면 calc_failed로 분류해서 숨기지 않음.
find_threshold_time_v2 <- function(times, concs, threshold) {
  if (length(times) < 2 || is.na(threshold) || threshold <= 0) {
    return(list(t = NA_real_, status = "calc_failed"))
  }
  idx <- which(concs <= threshold)[1]
  if (is.na(idx)) {
    return(list(t = NA_real_, status = "censored"))       # 관찰 종료까지 미도달 (우측 검열)
  }
  if (idx == 1) {
    return(list(t = NA_real_, status = "calc_failed"))     # 이론상 발생하면 안 되는 케이스 — 조용히 숨기지 않음
  }
  t1 <- times[idx - 1]; t2 <- times[idx]
  c1 <- concs[idx - 1]; c2 <- concs[idx]
  list(t = t1 + (threshold - c1) * (t2 - t1) / (c2 - c1), status = "reached")
}

trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# 환자 1명 PK metrics 추출 (v2)
# .x: group_modify 내부 data (post_time, conc, legacy_Css 컬럼 포함)
# C0: 중단 시점 순간 농도 (post_time=0에서의 conc와 동일해야 함)
extract_patient_metrics_v2 <- function(.x, C0, legacy_Css, t_half_impl_parent, t_half_impl_met, obs_h) {

  times <- .x$post_time
  concs <- .x$conc

  base_na <- tibble(
    C0 = NA_real_, legacy_Css = legacy_Css,
    reldiff_24h = NA_real_, reldiff_72h = NA_real_,
    pct_drop_24h = NA_real_, pct_drop_72h = NA_real_,
    t_75pct = NA_real_, status_75pct = "calc_failed",
    t_50pct = NA_real_, status_50pct = "calc_failed",
    t_25pct = NA_real_, status_25pct = "calc_failed",
    t_10pct = NA_real_, status_10pct = "calc_failed",
    AUC_deficit_7d = NA_real_, AUC_deficit_7d_norm_pct = NA_real_,
    t_half_impl_parent = t_half_impl_parent, t_half_impl_met = t_half_impl_met
  )

  if (is.na(C0) || C0 <= 0 || length(times) < 3) return(base_na)

  C_24h <- approx(times, concs, xout = 24)$y
  C_72h <- approx(times, concs, xout = 72)$y

  reldiff_24h <- if (!is.na(C_24h)) (C_24h - C0) / (24 * C0) else NA_real_
  reldiff_72h <- if (!is.na(C_72h)) (C_72h - C0) / (72 * C0) else NA_real_
  pct_drop_24h <- if (!is.na(C_24h)) (C0 - C_24h) / C0 * 100 else NA_real_
  pct_drop_72h <- if (!is.na(C_72h)) (C0 - C_72h) / C0 * 100 else NA_real_

  r75 <- find_threshold_time_v2(times, concs, 0.75 * C0)
  r50 <- find_threshold_time_v2(times, concs, 0.50 * C0)
  r25 <- find_threshold_time_v2(times, concs, 0.25 * C0)
  r10 <- find_threshold_time_v2(times, concs, 0.10 * C0)

  # AUC deficit — t50% 도달 여부와 무관하게 항상 계산 (7일 또는 obs_h 중 짧은 쪽)
  obs_max <- min(168, obs_h)
  obs_7d  <- min(obs_max, max(times, na.rm = TRUE))
  post_7d <- .x %>% filter(post_time <= obs_7d)
  if (nrow(post_7d) >= 2) {
    AUC_actual     <- trapz(post_7d$post_time, post_7d$conc)
    AUC_expected   <- C0 * obs_7d
    AUC_deficit_7d <- AUC_expected - AUC_actual
    AUC_deficit_7d_norm_pct <- 100 * AUC_deficit_7d / AUC_expected   # 환자별 정규화 (A2 단위 오류 수정)
  } else {
    AUC_deficit_7d <- NA_real_
    AUC_deficit_7d_norm_pct <- NA_real_
  }

  tibble(
    C0                       = round(C0, 6),
    legacy_Css               = round(legacy_Css, 6),
    reldiff_24h              = round(reldiff_24h, 6),
    reldiff_72h              = round(reldiff_72h, 6),
    pct_drop_24h             = round(pct_drop_24h, 3),
    pct_drop_72h             = round(pct_drop_72h, 3),
    t_75pct                  = round(r75$t, 3), status_75pct = r75$status,
    t_50pct                  = round(r50$t, 3), status_50pct = r50$status,
    t_25pct                  = round(r25$t, 3), status_25pct = r25$status,
    t_10pct                  = round(r10$t, 3), status_10pct = r10$status,
    AUC_deficit_7d           = round(AUC_deficit_7d, 4),
    AUC_deficit_7d_norm_pct  = round(AUC_deficit_7d_norm_pct, 4),
    t_half_impl_parent       = round(t_half_impl_parent, 3),
    t_half_impl_met          = if (!is.na(t_half_impl_met)) round(t_half_impl_met, 3) else NA_real_
  )
}

# ── 4. 배치 시뮬레이션 함수 (v2) ──────────────────────────────────────────────

run_batch_v2 <- function(drug_name, batch_df, cfg) {

  t_disc <- cfg$n_doses * cfg$tau_h

  ev <- et(amt  = cfg$F * cfg$dose_mg,   # v2: F=1이므로 전체 공칭 용량이 들어감
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

  # legacy_Css: v1과 동일한 방식(마지막 투약 간격 평균) — S1 비교용으로만 보존
  legacy_css_df <- sim %>%
    filter(time >= (cfg$n_doses - 1) * cfg$tau_h, time <= t_disc) %>%
    group_by(id) %>%
    summarise(legacy_Css = mean(conc, na.rm = TRUE), .groups = "drop")

  # C0: 중단 시점(t_disc)과 가장 가까운 격자점의 농도 — v2의 주 기준
  c0_df <- sim %>%
    group_by(id) %>%
    slice(which.min(abs(time - t_disc))) %>%
    ungroup() %>%
    transmute(id, C0 = conc)

  post_df <- sim %>%
    filter(time >= t_disc) %>%
    mutate(post_time = time - t_disc) %>%
    select(id, post_time, conc) %>%
    left_join(c0_df, by = "id") %>%
    left_join(legacy_css_df, by = "id")

  # 환자별 구현 반감기 (D2, Table 4 재분석용)
  half_life_df <- batch_df %>%
    transmute(
      id,
      t_half_impl_parent = 0.693 * Vd_sampled / CL_sampled,
      t_half_impl_met = if (cfg$model == "parent_met") 0.693 * Vd_met_sampled / CL_met_sampled else NA_real_
    )

  metrics <- post_df %>%
    left_join(half_life_df, by = "id") %>%
    group_by(id) %>%
    group_modify(~ extract_patient_metrics_v2(
      .x,
      C0 = unique(.x$C0),
      legacy_Css = unique(.x$legacy_Css),
      t_half_impl_parent = unique(.x$t_half_impl_parent),
      t_half_impl_met = unique(.x$t_half_impl_met),
      obs_h = cfg$obs_h
    )) %>%
    ungroup()

  batch_df %>%
    select(id, drug, cyp2d6, CL_sampled, Vd_sampled) %>%
    left_join(metrics, by = "id")
}

# ── 5. 전체 실행 ──────────────────────────────────────────────────────────────
# v1과의 핵심 차이: 새로 샘플링하지 않고 6-1단계에서 백업한 기존 코호트를
# 그대로 읽음 — 정의 변경(D1) 효과와 표본추출 차이를 분리하기 위함.

params_path <- file.path(BASE, "outputs", "mc", sprintf("mc_params_all_n%d_v1.csv", N))
if (!file.exists(params_path)) {
  stop(sprintf(
    "%s 없음. N=%d용 백업 파일이 없으면: (1) N=100으로 테스트할 땐 outputs/mc/mc_params_all_n100.csv를 _v1.csv로 복사해서 쓰거나, (2) N=10000 본 실행 전엔 6-1단계에서 만든 mc_params_all_n10000_v1.csv를 그대로 써야 함(새로 샘플링하면 안 됨).",
    params_path, N
  ))
}
params_all <- read_csv(params_path, show_col_types = FALSE)

results_list <- list()

for (drug_name in names(drug_config)) {

  cfg         <- drug_config[[drug_name]]
  drug_params <- params_all %>% filter(drug == drug_name)
  n_batches   <- ceiling(nrow(drug_params) / BATCH_SIZE)

  cat(sprintf("\n▶ [v2] %s: %d patients / %d batches (obs_h=%g h, dt=%g h)\n",
              drug_name, nrow(drug_params), n_batches, cfg$obs_h, cfg$dt))

  t_start      <- proc.time()["elapsed"]
  drug_results <- vector("list", n_batches)

  for (b in seq_len(n_batches)) {
    idx_start <- (b - 1) * BATCH_SIZE + 1
    idx_end   <- min(b * BATCH_SIZE, nrow(drug_params))
    batch_df  <- drug_params[idx_start:idx_end, ]
    drug_results[[b]] <- run_batch_v2(drug_name, batch_df, cfg)
    cat(sprintf("  batch %02d/%02d done (n=%d)\n", b, n_batches, nrow(batch_df)))
  }

  drug_df <- bind_rows(drug_results)
  elapsed <- round(proc.time()["elapsed"] - t_start, 1)

  out_path <- file.path(OUT_DIR, sprintf("mc_results_v2_%s_n%d.csv",
                                          tolower(gsub("_", "", drug_name)), N))
  write_csv(drug_df, out_path)
  cat(sprintf("✓ Saved: %s  [%.1f sec]\n", out_path, elapsed))

  results_list[[drug_name]] <- drug_df
}

# ── 6. 전체 합본 저장 ─────────────────────────────────────────────────────────

results_all <- bind_rows(results_list)
all_path    <- file.path(OUT_DIR, sprintf("mc_results_v2_all_n%d.csv", N))
write_csv(results_all, all_path)

cat(sprintf("\n✅ mc_run_v2.R 완료\n"))
cat(sprintf("   Total rows : %d\n", nrow(results_all)))
cat(sprintf("   Saved      : %s\n", all_path))

# ── 7. 빠른 검증 — status 코드별 분포 + 지표 요약 ─────────────────────────────

cat("\n── Status code distribution (t_50pct) ──\n")
results_all %>%
  count(drug, status_50pct) %>%
  group_by(drug) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  print(n = Inf)

cat("\n── Quick validation (median metrics by drug, reached 환자만) ──\n")
results_all %>%
  filter(status_50pct == "reached") %>%
  group_by(drug) %>%
  summarise(
    C0_ngmL         = round(median(C0, na.rm = TRUE) * 1000, 1),
    legacy_Css_ngmL = round(median(legacy_Css, na.rm = TRUE) * 1000, 1),
    t50_h           = round(median(t_50pct, na.rm = TRUE), 1),
    n_reached       = n(),
    AUCdef_norm_pct = round(median(AUC_deficit_7d_norm_pct, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  print()

cat("\n▶ 다음 단계: qualification_v2.R, sensitivity_v2.R, 표 생성 스크립트\n")
cat("   (N=100으로 에러 없이 끝까지 돌았으면 N을 10000으로 바꿔서 본 실행)\n")
