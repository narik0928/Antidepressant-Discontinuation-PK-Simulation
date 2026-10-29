# ============================================================
# sensitivity_v3.R -- Sensitivity Analysis v3
# Author: Nari Kim (code assist: Claude)
#
# HANDOFF_TO_NARI_v3.md 7-1-4 구현.
#
# v2 대비 바뀐 것
#   fm,CYP2D6 sweep (항목 3) -- 없앰(D5: fluoxetine.venlafaxine은 더 이상
#     CYP2D6를 모델링하지 않으므로 이 민감도 자체가 의미 없어짐).
#   신규: paroxetine PM 배수 민감도(H14) -- 0.58(주분석, Sindrup 1992)과
#     0.40(Liao 2024 보정 Css비 2.50의 역수) 두 조건을 비교. 조건별로
#     (1) PM 군의 중앙값 t50%, (2) censoring 수, (3) PM 대 non-PM
#     효과크기(rank-biserial r, Mann-Whitney U 기반)를 출력.
#   potency weight sweep (항목 5) -- fluoxetine, venlafaxine(IR/XR)에서
#     유지. v3는 CYP2D6를 이 약물들에 적용하지 않으므로(D5) 기반 샘플링에서
#     phenotype 배정 자체를 제거(v2의 apply_fm_cyp2d6 불필요). D11
#     중앙값보존 로그정규 샘플링을 씀.
#   파라미터: pk_parameters_v3.R에서 읽음(하드코딩 금지 원칙 유지).
#   추가 지시 1(2-1) 반영: paroxetine V/F CV 민감도(0.35 주분석 vs 0.835
#     Li 2022) 신규 추가 -- t50%/t25%/t10% censoring, PM/non-PM 비·효과크기.
#   추가 지시 1(2-2) 반영: venlafaxine IR/XR 관찰창 7일->14일(drug_config).
#
# 사용법: 먼저 N_sens=500으로 에러 없이 끝까지 도는지 확인 -> 필요시
#   2000으로 올려 부트스트랩 CI 폭이 좁아지는지 재확인.
# ============================================================

suppressPackageStartupMessages({
  library(rxode2)
  library(tidyverse)
})

BASE    <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR <- file.path(BASE, "outputs", "mc_sensitivity_v3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "data/pk_parameters/pk_parameters_v3.R"))   # v3: D2'/D5/D7/D9/D11 반영본

set.seed(2025)

N_sens <- 2000   # 먼저 500으로 테스트

# ------------------------------------------------------------
# 1. rxode2 모델 (v2와 동일 구조)
# ------------------------------------------------------------
ode_1cmt <- rxode2({
  Cp = A_c / Vd
  d/dt(depot) = -Ka * depot
  d/dt(A_c)   =  Ka * depot - (CL / Vd) * A_c
})

ode_parent_met <- rxode2({
  Cp = A_p / Vd_p
  Cm = A_m / Vd_m
  d/dt(depot) = -Ka * depot
  d/dt(A_p)   =  Ka * depot - (CL_p / Vd_p) * A_p
  d/dt(A_m)   =  Fm * CL_p * Cp - (CL_m / Vd_m) * A_m
})

# ------------------------------------------------------------
# 2. drug_config (v3, HANDOFF_v3 6절 투여설계)
# ------------------------------------------------------------
drug_config <- list(
  fluoxetine = list(
    model = "parent_met",
    Ka = pk_params$fluoxetine$Ka,
    F  = pk_params$fluoxetine$F,
    Fm = pk_params$fluoxetine$norfluoxetine$Fm,
    dose_mg = 20, tau_h = 24, n_doses = 60, obs_h = 400 * 24, dt = 2.0
  ),
  venlafaxine_IR = list(
    model = "parent_met",
    Ka = pk_params$venlafaxine_IR$Ka,
    F  = pk_params$venlafaxine_IR$F,
    Fm = pk_params$venlafaxine_IR$ODV$Fm,
    dose_mg = 25, tau_h = 8, n_doses = 21, obs_h = 14 * 24, dt = 0.25  # 추가 지시 1(2-2): 7일->14일
  ),
  venlafaxine_XR = list(
    model = "parent_met",
    Ka = pk_params$venlafaxine_XR$Ka,
    F  = pk_params$venlafaxine_XR$F,
    Fm = pk_params$venlafaxine_XR$ODV$Fm,
    dose_mg = 75, tau_h = 24, n_doses = 7, obs_h = 14 * 24, dt = 0.25  # 추가 지시 1(2-2): 7일->14일
  )
)

paroxetine_cfg <- list(
  Ka = pk_params$paroxetine$Ka,
  F  = pk_params$paroxetine$F,
  dose_mg = 30, tau_h = 24, n_doses = 14, obs_h = 30 * 24, dt = 0.5  # D6/B6
)

# ------------------------------------------------------------
# 3. D11 중앙값보존 로그정규 샘플링 (v3)
# ------------------------------------------------------------
sample_lognormal_v3 <- function(n, pop_value, cv) {
  sigma <- sqrt(log(cv^2 + 1))
  mu    <- log(pop_value)              # D11: 중앙값보존 (- sigma^2/2 없음)
  rlnorm(n, meanlog = mu, sdlog = sigma)
}

# 항목 5(potency) 기반 샘플링 -- CYP2D6 배정 없음(D5: 이 세 약물은 모델링 안 함)
sample_drug_base_v3 <- function(n, drug_key) {
  p <- pk_params[[drug_key]]
  base <- tibble(
    id         = seq_len(n),
    CL_sampled = sample_lognormal_v3(n, p$CL, p$CV_CL),
    Vd_sampled = sample_lognormal_v3(n, p$Vd, p$CV_Vd)
  )
  if (drug_key == "fluoxetine") {
    met <- p$norfluoxetine
    base <- base %>% mutate(
      CL_met_sampled = sample_lognormal_v3(n, met$CL_met, met$CV_CL_met),
      Vd_met_sampled = rep(met$Vd_met, n)
    )
  }
  if (drug_key %in% c("venlafaxine_IR", "venlafaxine_XR")) {
    met <- p$ODV
    base <- base %>% mutate(
      CL_met_sampled = sample_lognormal_v3(n, met$CL_odv, met$CV_CL_odv),
      Vd_met_sampled = rep(met$Vd_odv, n)
    )
  }
  base
}

# ------------------------------------------------------------
# 4. 상태 코드가 있는 역치-도달 시각 함수 (v2와 동일 로직)
# ------------------------------------------------------------
find_threshold_time_v3 <- function(times, concs, threshold) {
  if (length(times) < 2 || is.na(threshold) || threshold <= 0) {
    return(list(t = NA_real_, status = "calc_failed"))
  }
  idx <- which(concs <= threshold)[1]
  if (is.na(idx)) return(list(t = NA_real_, status = "censored"))
  if (idx == 1)    return(list(t = NA_real_, status = "calc_failed"))
  t1 <- times[idx - 1]; t2 <- times[idx]
  c1 <- concs[idx - 1]; c2 <- concs[idx]
  list(t = t1 + (threshold - c1) * (t2 - t1) / (c2 - c1), status = "reached")
}

trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

bootstrap_median_ci <- function(x, B = 2000, conf = 0.95) {
  x <- x[!is.na(x)]
  if (length(x) < 5) return(c(lo = NA_real_, hi = NA_real_))
  boots <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  alpha <- (1 - conf) / 2
  c(lo = unname(quantile(boots, alpha)), hi = unname(quantile(boots, 1 - alpha)))
}

# rank-biserial r from Mann-Whitney U (two independent groups)
rank_biserial_r <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  w <- suppressWarnings(wilcox.test(x1, x2, exact = FALSE))
  U <- w$statistic
  n1 <- length(x1); n2 <- length(x2)
  1 - (2 * U) / (n1 * n2)
}

# ============================================================
# 5. 항목 5 -- potency weight sweep (fluoxetine, venlafaxine IR/XR)
# ============================================================
potency_grid   <- c(0.5, 0.75, 1.0, 1.25, 1.5)
drugs_potency  <- c("fluoxetine", "venlafaxine_IR", "venlafaxine_XR")

run_potency_batch_v3 <- function(params_df, cfg, potency_weight) {
  t_disc <- cfg$n_doses * cfg$tau_h
  ev <- et(amt = cfg$F * cfg$dose_mg, ii = cfg$tau_h, addl = cfg$n_doses - 1, time = 0) %>%
        et(seq(0, t_disc + cfg$obs_h, by = cfg$dt))

  inits <- c(depot = 0, A_p = 0, A_m = 0)
  get_p <- function(row) c(Ka = cfg$Ka, Vd_p = row$Vd_sampled, CL_p = row$CL_sampled,
                            Fm = cfg$Fm, Vd_m = row$Vd_met_sampled, CL_m = row$CL_met_sampled)

  sim_list <- lapply(seq_len(nrow(params_df)), function(i) {
    row <- params_df[i, , drop = FALSE]
    sim_i <- rxSolve(ode_parent_met, get_p(row), ev, inits = inits) %>% as.data.frame()
    sim_i$id <- row$id
    sim_i
  })
  sim <- bind_rows(sim_list) %>% mutate(conc = Cp + potency_weight * Cm)

  c0_df <- sim %>% group_by(id) %>% slice(which.min(abs(time - t_disc))) %>%
    ungroup() %>% transmute(id, C0 = conc)

  post_df <- sim %>% filter(time >= t_disc) %>%
    mutate(post_time = time - t_disc) %>%
    select(id, post_time, conc) %>%
    left_join(c0_df, by = "id")

  obs_7d_cap <- min(168, cfg$obs_h)

  post_df %>%
    group_by(id) %>%
    group_modify(~{
      times <- .x$post_time; concs <- .x$conc; C0 <- unique(.x$C0)
      r50 <- find_threshold_time_v3(times, concs, 0.5 * C0)

      obs_7d  <- min(obs_7d_cap, max(times, na.rm = TRUE))
      post_7d <- .x %>% filter(post_time <= obs_7d)
      if (nrow(post_7d) >= 2 && !is.na(C0) && C0 > 0) {
        auc_actual   <- trapz(post_7d$post_time, post_7d$conc)
        auc_expected <- C0 * obs_7d
        auc_def_abs  <- auc_expected - auc_actual
        auc_def_norm <- 100 * auc_def_abs / auc_expected
      } else {
        auc_def_abs  <- NA_real_; auc_def_norm <- NA_real_
      }

      tibble(C0 = C0, t50_C0 = r50$t, status_t50 = r50$status,
             AUC_deficit_7d_abs = auc_def_abs, AUC_deficit_7d_norm_pct = auc_def_norm)
    }) %>%
    ungroup()
}

item5_results <- list()
for (drug in drugs_potency) {
  cfg <- drug_config[[drug]]
  base_df <- sample_drug_base_v3(N_sens, drug)
  for (w in potency_grid) {
    cat(sprintf("[v3 potency] %s | potency_weight=%.2f ...\n", drug, w))
    res <- run_potency_batch_v3(base_df, cfg, potency_weight = w)
    res$drug <- drug
    res$potency_weight <- w
    item5_results[[paste(drug, w)]] <- res
  }
}
item5_all <- bind_rows(item5_results)

item5_summary <- item5_all %>%
  group_by(drug, potency_weight) %>%
  summarise(
    n_total          = n(),
    n_reached        = sum(status_t50 == "reached", na.rm = TRUE),
    n_censored       = sum(status_t50 == "censored", na.rm = TRUE),
    n_calc_failed    = sum(status_t50 == "calc_failed", na.rm = TRUE),
    median_t50_h     = median(t50_C0[status_t50 == "reached"], na.rm = TRUE),
    t50_ci_lo        = bootstrap_median_ci(t50_C0[status_t50 == "reached"])["lo"],
    t50_ci_hi        = bootstrap_median_ci(t50_C0[status_t50 == "reached"])["hi"],
    median_AUCdef_norm_pct = median(AUC_deficit_7d_norm_pct, na.rm = TRUE),
    AUCdef_norm_ci_lo = bootstrap_median_ci(AUC_deficit_7d_norm_pct)["lo"],
    AUCdef_norm_ci_hi = bootstrap_median_ci(AUC_deficit_7d_norm_pct)["hi"],
    .groups = "drop"
  ) %>%
  arrange(drug, potency_weight)

cat("\n========== v3 항목5: potency weight sweep ==========\n")
print(item5_summary, n = 100)
write_csv(item5_all,     file.path(OUT_DIR, "item5_potency_raw_v3.csv"))
write_csv(item5_summary, file.path(OUT_DIR, "item5_potency_summary_v3.csv"))

# ============================================================
# 6. 신규 -- paroxetine PM 배수 민감도 (H14)
#    0.58(주분석, Sindrup 1992) vs 0.40(Liao 2024 보정 Css비 2.50의 역수)
# ============================================================
pm_factor_grid <- c(main_0.58 = 0.58, alt_0.40 = 0.40)

run_paroxetine_batch_v3 <- function(params_df, cfg) {
  t_disc <- cfg$n_doses * cfg$tau_h
  ev <- et(amt = cfg$F * cfg$dose_mg, ii = cfg$tau_h, addl = cfg$n_doses - 1, time = 0) %>%
        et(seq(0, t_disc + cfg$obs_h, by = cfg$dt))
  inits <- c(depot = 0, A_c = 0)
  get_p <- function(row) c(Ka = cfg$Ka, Vd = row$Vd_sampled, CL = row$CL_sampled)

  sim_list <- lapply(seq_len(nrow(params_df)), function(i) {
    row <- params_df[i, , drop = FALSE]
    sim_i <- rxSolve(ode_1cmt, get_p(row), ev, inits = inits) %>% as.data.frame()
    sim_i$id <- row$id
    sim_i
  })
  sim <- bind_rows(sim_list) %>% mutate(conc = Cp)

  c0_df <- sim %>% group_by(id) %>% slice(which.min(abs(time - t_disc))) %>%
    ungroup() %>% transmute(id, C0 = conc)

  post_df <- sim %>% filter(time >= t_disc) %>%
    mutate(post_time = time - t_disc) %>%
    select(id, post_time, conc) %>%
    left_join(c0_df, by = "id") %>%
    left_join(params_df %>% select(id, cyp2d6), by = "id")

  post_df %>%
    group_by(id, cyp2d6) %>%
    group_modify(~{
      times <- .x$post_time; concs <- .x$conc; C0 <- unique(.x$C0)
      r50 <- find_threshold_time_v3(times, concs, 0.5 * C0)
      tibble(C0 = C0, t50_C0 = r50$t, status_t50 = r50$status)
    }) %>%
    ungroup()
}

item_pm_results <- list()
for (cond_name in names(pm_factor_grid)) {
  fac <- pm_factor_grid[[cond_name]]
  cat(sprintf("[v3 PM sensitivity] paroxetine | PM factor=%.2f (%s) ...\n", fac, cond_name))

  cyp <- sample(c("poor", "normal"), N_sens, replace = TRUE,
                prob = c(cyp2d6$poor_metabolizer$freq, cyp2d6$normal_metabolizer$freq))
  fac_vec <- if_else(cyp == "poor", fac, 1.0)

  params_df <- tibble(
    id         = seq_len(N_sens),
    cyp2d6     = cyp,
    CL_sampled = sample_lognormal_v3(N_sens, pk_params$paroxetine$CL, pk_params$paroxetine$CV_CL) * fac_vec,
    Vd_sampled = sample_lognormal_v3(N_sens, pk_params$paroxetine$Vd, pk_params$paroxetine$CV_Vd)
  )

  res <- run_paroxetine_batch_v3(params_df, paroxetine_cfg)
  res$pm_condition <- cond_name
  res$pm_factor <- fac
  item_pm_results[[cond_name]] <- res
}
item_pm_all <- bind_rows(item_pm_results)

item_pm_summary <- item_pm_all %>%
  group_by(pm_condition, pm_factor, cyp2d6) %>%
  summarise(
    n_total        = n(),
    n_reached      = sum(status_t50 == "reached", na.rm = TRUE),
    n_censored     = sum(status_t50 == "censored", na.rm = TRUE),
    pct_censored   = round(100 * n_censored / n_total, 2),
    median_t50_h   = median(t50_C0[status_t50 == "reached"], na.rm = TRUE),
    .groups = "drop"
  )

item_pm_effect <- item_pm_all %>%
  group_by(pm_condition, pm_factor) %>%
  summarise(
    effect_size_rank_biserial_r = rank_biserial_r(
      t50_C0[cyp2d6 == "poor"   & status_t50 == "reached"],
      t50_C0[cyp2d6 == "normal" & status_t50 == "reached"]
    ),
    wilcoxon_p = suppressWarnings(wilcox.test(
      t50_C0[cyp2d6 == "poor"   & status_t50 == "reached"],
      t50_C0[cyp2d6 == "normal" & status_t50 == "reached"],
      exact = FALSE
    )$p.value),
    .groups = "drop"
  )

cat("\n========== v3 신규: paroxetine PM 배수 민감도 (H14) ==========\n")
cat("-- group summary (PM median t50%, censoring) --\n")
print(item_pm_summary, n = Inf)
cat("\n-- PM vs non-PM effect size (rank-biserial r, Mann-Whitney 기반) --\n")
print(item_pm_effect, n = Inf)

write_csv(item_pm_all,     file.path(OUT_DIR, "item_pm_sensitivity_raw_v3.csv"))
write_csv(item_pm_summary, file.path(OUT_DIR, "item_pm_sensitivity_summary_v3.csv"))
write_csv(item_pm_effect,  file.path(OUT_DIR, "item_pm_sensitivity_effectsize_v3.csv"))

# ============================================================
# 6b. 신규(추가 지시 1 2-1) -- paroxetine V/F CV 민감도
#     0.35(주분석, 저자 가정) vs 0.835(Li 2022 Table 2)
#     조건별 t50%/t25%/t10% 중앙값·사분위수, PM 대 non-PM 중앙값 비,
#     PM 대 non-PM 효과크기(rank-biserial r), 지표별 censoring 수를
#     출력. 덕현 쪽 기대값(expected_mc_v3_parV835.csv/summary, 3-4
#     게이트 7번)과 대조.
# ============================================================
par_vd_cv_grid <- c(main_0.35 = 0.35, alt_0.835 = 0.835)

run_paroxetine_batch_v3_multi <- function(params_df, cfg) {
  t_disc <- cfg$n_doses * cfg$tau_h
  ev <- et(amt = cfg$F * cfg$dose_mg, ii = cfg$tau_h, addl = cfg$n_doses - 1, time = 0) %>%
        et(seq(0, t_disc + cfg$obs_h, by = cfg$dt))
  inits <- c(depot = 0, A_c = 0)
  get_p <- function(row) c(Ka = cfg$Ka, Vd = row$Vd_sampled, CL = row$CL_sampled)

  sim_list <- lapply(seq_len(nrow(params_df)), function(i) {
    row <- params_df[i, , drop = FALSE]
    sim_i <- rxSolve(ode_1cmt, get_p(row), ev, inits = inits) %>% as.data.frame()
    sim_i$id <- row$id
    sim_i
  })
  sim <- bind_rows(sim_list) %>% mutate(conc = Cp)

  c0_df <- sim %>% group_by(id) %>% slice(which.min(abs(time - t_disc))) %>%
    ungroup() %>% transmute(id, C0 = conc)

  post_df <- sim %>% filter(time >= t_disc) %>%
    mutate(post_time = time - t_disc) %>%
    select(id, post_time, conc) %>%
    left_join(c0_df, by = "id") %>%
    left_join(params_df %>% select(id, cyp2d6), by = "id")

  post_df %>%
    group_by(id, cyp2d6) %>%
    group_modify(~{
      times <- .x$post_time; concs <- .x$conc; C0 <- unique(.x$C0)
      r50 <- find_threshold_time_v3(times, concs, 0.50 * C0)
      r25 <- find_threshold_time_v3(times, concs, 0.25 * C0)
      r10 <- find_threshold_time_v3(times, concs, 0.10 * C0)
      tibble(C0 = C0,
             t50_C0 = r50$t, status_t50 = r50$status,
             t25_C0 = r25$t, status_t25 = r25$status,
             t10_C0 = r10$t, status_t10 = r10$status)
    }) %>%
    ungroup()
}

item_parvdcv_results <- list()
for (cond_name in names(par_vd_cv_grid)) {
  cv_v <- par_vd_cv_grid[[cond_name]]
  cat(sprintf("[v3 V/F CV sensitivity] paroxetine | CV_Vd=%.3f (%s) ...\n", cv_v, cond_name))

  cyp <- sample(c("poor", "normal"), N_sens, replace = TRUE,
                prob = c(cyp2d6$poor_metabolizer$freq, cyp2d6$normal_metabolizer$freq))
  fac_vec <- if_else(cyp == "poor", cyp2d6$poor_metabolizer$CL_factor, 1.0)

  params_df <- tibble(
    id         = seq_len(N_sens),
    cyp2d6     = cyp,
    CL_sampled = sample_lognormal_v3(N_sens, pk_params$paroxetine$CL, pk_params$paroxetine$CV_CL) * fac_vec,
    Vd_sampled = sample_lognormal_v3(N_sens, pk_params$paroxetine$Vd, cv_v)
  )

  res <- run_paroxetine_batch_v3_multi(params_df, paroxetine_cfg)
  res$vdcv_condition <- cond_name
  res$vdcv_value <- cv_v
  item_parvdcv_results[[cond_name]] <- res
}
item_parvdcv_all <- bind_rows(item_parvdcv_results)

item_parvdcv_summary <- item_parvdcv_all %>%
  group_by(vdcv_condition, vdcv_value) %>%
  summarise(
    n_total           = n(),
    median_t50_h      = median(t50_C0[status_t50 == "reached"], na.rm = TRUE),
    q1_t50_h          = quantile(t50_C0[status_t50 == "reached"], 0.25, na.rm = TRUE),
    q3_t50_h          = quantile(t50_C0[status_t50 == "reached"], 0.75, na.rm = TRUE),
    n_censored_t50    = sum(status_t50 == "censored", na.rm = TRUE),
    n_censored_t25    = sum(status_t25 == "censored", na.rm = TRUE),
    n_censored_t10    = sum(status_t10 == "censored", na.rm = TRUE),
    median_t50_PM     = median(t50_C0[cyp2d6 == "poor"   & status_t50 == "reached"], na.rm = TRUE),
    median_t50_nonPM  = median(t50_C0[cyp2d6 == "normal" & status_t50 == "reached"], na.rm = TRUE),
    ratio_PM_nonPM    = median_t50_PM / median_t50_nonPM,
    .groups = "drop"
  )

item_parvdcv_effect <- item_parvdcv_all %>%
  group_by(vdcv_condition, vdcv_value) %>%
  summarise(
    effect_size_rank_biserial_r = rank_biserial_r(
      t50_C0[cyp2d6 == "poor"   & status_t50 == "reached"],
      t50_C0[cyp2d6 == "normal" & status_t50 == "reached"]
    ),
    wilcoxon_p = suppressWarnings(wilcox.test(
      t50_C0[cyp2d6 == "poor"   & status_t50 == "reached"],
      t50_C0[cyp2d6 == "normal" & status_t50 == "reached"],
      exact = FALSE
    )$p.value),
    .groups = "drop"
  )

cat("\n========== v3 신규(추가 지시 1 2-1): paroxetine V/F CV 민감도 ==========\n")
cat("-- summary (t50% 중앙값/사분위수, PM/non-PM 비, censoring) --\n")
print(item_parvdcv_summary, n = Inf)
cat("\n-- PM vs non-PM effect size (rank-biserial r) --\n")
print(item_parvdcv_effect, n = Inf)
cat("\n기대값(덕현 쪽 독립 계산, HANDOFF_v3_ADDENDUM_1028 2-1): 0.35 조건 t50% 중앙값 21.9h\n")
cat("[14.4-33.2], PM/non-PM비 1.64, r=0.458; 0.835 조건 t50% 중앙값 22.0h [12.1-39.8],\n")
cat("PM/non-PM비 1.66, r=0.327. censoring(0.35/0.835): t50% 0/2, t25% 0/9, t10% 2/52.\n")
cat("self_check_v3.py --expected-mc expected_mc_v3_parV835.csv --expected-params\n")
cat("expected_params_v3_parV835.csv 로 N=10,000 본 실행 결과와 대조할 것 (3-4 게이트 7번).\n")

write_csv(item_parvdcv_all,     file.path(OUT_DIR, "item_par_vdcv_sensitivity_raw_v3.csv"))
write_csv(item_parvdcv_summary, file.path(OUT_DIR, "item_par_vdcv_sensitivity_summary_v3.csv"))
write_csv(item_parvdcv_effect,  file.path(OUT_DIR, "item_par_vdcv_sensitivity_effectsize_v3.csv"))

# ============================================================
# 7. N 적정성 메모
# ============================================================
cat("\n-- N_sens 적정성 점검 --\n")
cat(sprintf("현재 N_sens = %d. ci_lo~ci_hi 폭이 median 대비 넓으면 N_sens를 2000(이상)으로\n", N_sens))
cat("올려서 다시 실행.\n")

cat("\n완료. outputs/mc_sensitivity_v3/ 아래 8개 CSV 생성됨.\n")
cat("(fm,CYP2D6 sweep은 D5에 따라 v3에서 제거됨 -- 대신 paroxetine PM 배수 민감도 +\n")
cat(" paroxetine V/F CV 민감도[추가 지시 1 2-1]로 대체)\n\n")
