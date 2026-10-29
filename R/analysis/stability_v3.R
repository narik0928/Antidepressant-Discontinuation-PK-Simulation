# ============================================================
# stability_v3.R -- N/seed Stability Check + Fluoxetine Interpolation
# Error Check (v3)
# Author: Nari Kim (code assist: Claude)
#
# HANDOFF_TO_NARI_v3.md 7-1-5 구현 (A11, H16):
#   "주 분석 중앙값과 효과크기가 N과 seed에 안정적인지 확인합니다. 예:
#   N = 2,000, 5,000, 10,000과 seed 3개. fluoxetine 2h 시간격자의 보간
#   오차도 몇 명 표본에서 0.5h 격자와 비교하세요."
#
# 설계 메모: 기존 R/analysis/stability_check.R(v1)은 이미 뽑힌 n=10000
# 결과에서 "서브샘플링"만 했음(재시뮬 없음). HANDOFF v3가 요구하는 "N과
# seed에 안정적인지"는 서로 다른 seed로 독립적으로 다시 뽑고 다시 돌려야
# 확인되는 것이라, 이 스크립트는 매 (N, seed) 조합마다 독립적으로 새로
# 샘플링(mc_params_v3.R의 fresh-sample 방식과 동일한 D11 중앙값보존
# 로그정규 + paroxetine만 PM x0.58)하고 ODE를 다시 푼다. 본 분석(z값 역산,
# N=10000, 결정론적)과는 별개의 확인 절차 -- "그 N=10000 결과가 샘플링
# 변동에 민감하지 않다"는 것을 보이기 위함.
#
# 핵심 지표만 계산(stability_check.R처럼 가볍게): 약물별 median t50%
# (reached-only), t25% censoring 비율, paroxetine PM 대 non-PM 효과크기
# (rank-biserial r, sensitivity_v3.R과 동일 지표).
#
# 주의: N=10000 x seed 3개 x drug 5개를 전부 새로 ODE로 도는 것이라
# 시간이 걸림(본 실행과 맞먹는 양). 먼저 N=2000만으로 끝까지 도는지
# 확인한 뒤 전체 실행 권장.
#
# Output: outputs/stability_v3/stability_cv_v3.csv
#         outputs/stability_v3/fluoxetine_interp_error_v3.csv
# ============================================================

suppressPackageStartupMessages({
  library(rxode2)
  library(tidyverse)
})

BASE    <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR <- file.path(BASE, "outputs", "stability_v3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "data/pk_parameters/pk_parameters_v3.R"))

N_GRID    <- c(2000, 5000, 10000)   # 먼저 2000만으로 테스트 권장
SEED_GRID <- c(2025, 2026, 2027)    # 주 분석 seed(2025) + 추가 2개

# ── 1. ODE 모델 ───────────────────────────────────────────────────────────────
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

# ── 2. drug_config (v3, mc_run_v3.R과 동일) ──────────────────────────────────
# 2026-10-28 수정: obs_h를 추가 지시 1(2-2)의 mc_run_v3.R/sensitivity_v3.R
# 보정값과 맞춤 -- 이 스크립트는 그동안 addendum 이전(구) 창(sertraline
# 14일, venlafaxine 7일)으로 남아있었음(실행 중 발견, 본 분석과 어긋나던
# 버그). Fluoxetine/Paroxetine은 원래도 addendum 영향 없음/이미 일치.
drug_config <- list(
  Sertraline = list(model = "1cmt", Ka = pk_params$sertraline$Ka, F = pk_params$sertraline$F,
                     dose_mg = 50, tau_h = 24, n_doses = 14, obs_h = 30 * 24, dt = 1.0, conc_col = "Cp"),
  Fluoxetine = list(model = "parent_met", Ka = pk_params$fluoxetine$Ka, F = pk_params$fluoxetine$F,
                     Fm = pk_params$fluoxetine$norfluoxetine$Fm,
                     dose_mg = 20, tau_h = 24, n_doses = 60, obs_h = 400 * 24, dt = 2.0, conc_col = "Ctotal"),
  Paroxetine = list(model = "1cmt", Ka = pk_params$paroxetine$Ka, F = pk_params$paroxetine$F,
                     dose_mg = 30, tau_h = 24, n_doses = 14, obs_h = 30 * 24, dt = 0.5, conc_col = "Cp"),
  Venlafaxine_IR = list(model = "parent_met", Ka = pk_params$venlafaxine_IR$Ka, F = pk_params$venlafaxine_IR$F,
                     Fm = pk_params$venlafaxine_IR$ODV$Fm,
                     dose_mg = 25, tau_h = 8, n_doses = 21, obs_h = 14 * 24, dt = 0.25, conc_col = "Ctotal"),
  Venlafaxine_XR = list(model = "parent_met", Ka = pk_params$venlafaxine_XR$Ka, F = pk_params$venlafaxine_XR$F,
                     Fm = pk_params$venlafaxine_XR$ODV$Fm,
                     dose_mg = 75, tau_h = 24, n_doses = 7, obs_h = 14 * 24, dt = 0.25, conc_col = "Ctotal")
)

# 환자 1명당 ODE 풀이 제한시간(초) -- 극단적으로 샘플링된 CL/Vd 때문에
# 솔버가 멈추는 경우를 막기 위함(2026-10-28, N=5000/seed=2026/Fluoxetine에서
# 1시간 넘게 멎는 현상 발견 후 추가). 정상 환자는 수십~수백 ms면 끝남 --
# 10초 넘으면 사실상 비정상으로 간주. 타임아웃/에러가 난 환자는 건너뛰고
# PATHOLOGICAL_LOG에 파라미터를 남긴다(결과에서는 제외, 조용히 버리지 않음).
PATIENT_TIMEOUT_SEC <- 10
PATHOLOGICAL_LOG <- list()

safe_rxSolve_v3 <- function(model, params, ev, inits, meta) {
  out <- tryCatch({
    setTimeLimit(elapsed = PATIENT_TIMEOUT_SEC, transient = TRUE)
    res <- as.data.frame(rxSolve(model, params, ev, inits = inits))
    setTimeLimit(elapsed = Inf, transient = TRUE)
    res
  }, error = function(e) {
    setTimeLimit(elapsed = Inf, transient = TRUE)
    PATHOLOGICAL_LOG[[length(PATHOLOGICAL_LOG) + 1]] <<- c(meta, list(error = conditionMessage(e)))
    cat(sprintf("  [WARN] 환자 건너뜀 (%s, id=%s): %s | params=%s\n",
                meta$drug, meta$id, conditionMessage(e),
                paste(names(params), round(params, 4), sep = "=", collapse = ", ")))
    NULL
  })
  out
}

sample_lognormal_v3 <- function(n, pop_value, cv) {
  sigma <- sqrt(log(cv^2 + 1)); mu <- log(pop_value)   # D11
  rlnorm(n, meanlog = mu, sdlog = sigma)
}

find_threshold_time_v3 <- function(times, concs, threshold) {
  if (length(times) < 2 || is.na(threshold) || threshold <= 0) return(list(t = NA_real_, status = "calc_failed"))
  idx <- which(concs <= threshold)[1]
  if (is.na(idx)) return(list(t = NA_real_, status = "censored"))
  if (idx == 1)    return(list(t = NA_real_, status = "calc_failed"))
  t1 <- times[idx - 1]; t2 <- times[idx]; c1 <- concs[idx - 1]; c2 <- concs[idx]
  list(t = t1 + (threshold - c1) * (t2 - t1) / (c2 - c1), status = "reached")
}

rank_biserial_r <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  U <- suppressWarnings(wilcox.test(x1, x2, exact = FALSE))$statistic
  1 - (2 * U) / (length(x1) * length(x2))
}

# 약물 1개, N명 샘플링 + ODE 실행 + t50%/t25% 추출 (가벼운 버전)
sample_and_run <- function(drug_name, n, seed) {
  set.seed(seed)
  cfg <- drug_config[[drug_name]]

  if (drug_name == "Paroxetine") {
    cyp <- sample(c("poor", "normal"), n, replace = TRUE,
                  prob = c(cyp2d6$poor_metabolizer$freq, cyp2d6$normal_metabolizer$freq))
    fac <- if_else(cyp == "poor", cyp2d6$poor_metabolizer$CL_factor, 1.0)
    params_df <- tibble(
      id = seq_len(n), cyp2d6 = cyp,
      CL_sampled = sample_lognormal_v3(n, pk_params$paroxetine$CL, pk_params$paroxetine$CV_CL) * fac,
      Vd_sampled = sample_lognormal_v3(n, pk_params$paroxetine$Vd, pk_params$paroxetine$CV_Vd)
    )
  } else {
    drug_key <- tolower(gsub("_IR|_XR", "", drug_name))
    drug_key <- if (grepl("_IR", drug_name)) "venlafaxine_IR" else if (grepl("_XR", drug_name)) "venlafaxine_XR" else drug_key
    p <- pk_params[[drug_key]]
    params_df <- tibble(
      id = seq_len(n), cyp2d6 = "not_modeled",
      CL_sampled = sample_lognormal_v3(n, p$CL, p$CV_CL),
      Vd_sampled = sample_lognormal_v3(n, p$Vd, p$CV_Vd)
    )
    if (cfg$model == "parent_met") {
      met_CL <- if (!is.null(p$ODV)) p$ODV$CL_odv else p$norfluoxetine$CL_met
      met_CV <- if (!is.null(p$ODV)) p$ODV$CV_CL_odv else p$norfluoxetine$CV_CL_met
      met_Vd <- if (!is.null(p$ODV)) p$ODV$Vd_odv else p$norfluoxetine$Vd_met
      params_df <- params_df %>% mutate(
        CL_met_sampled = sample_lognormal_v3(n, met_CL, met_CV),
        Vd_met_sampled = met_Vd
      )
    }
  }

  # 진단: 샘플링 직후 CL/Vd 범위 출력 (극단값을 ODE 풀이 전에 바로 확인)
  cat(sprintf("    CL_sampled range: [%.4g, %.4g] | Vd_sampled range: [%.4g, %.4g]\n",
              min(params_df$CL_sampled), max(params_df$CL_sampled),
              min(params_df$Vd_sampled), max(params_df$Vd_sampled)))
  if ("CL_met_sampled" %in% names(params_df)) {
    cat(sprintf("    CL_met_sampled range: [%.4g, %.4g]\n",
                min(params_df$CL_met_sampled), max(params_df$CL_met_sampled)))
  }

  t_disc <- cfg$n_doses * cfg$tau_h
  ev <- et(amt = cfg$F * cfg$dose_mg, ii = cfg$tau_h, addl = cfg$n_doses - 1, time = 0) %>%
        et(seq(0, t_disc + cfg$obs_h, by = cfg$dt))

  if (cfg$model == "1cmt") {
    model <- ode_1cmt; inits <- c(depot = 0, Cp = 0)
    get_p <- function(row) c(Ka = cfg$Ka, Vd = row$Vd_sampled, CL = row$CL_sampled)
  } else {
    model <- ode_parent_met; inits <- c(depot = 0, Cp = 0, Cm = 0)
    get_p <- function(row) c(Ka = cfg$Ka, Vd_p = row$Vd_sampled, CL_p = row$CL_sampled,
                              Fm = cfg$Fm, Vd_m = row$Vd_met_sampled, CL_m = row$CL_met_sampled)
  }

  t_batch_start <- proc.time()["elapsed"]
  sim_list <- lapply(seq_len(nrow(params_df)), function(i) {
    row <- params_df[i, , drop = FALSE]
    if (i %% 500 == 0) {
      cat(sprintf("    ... %d/%d patients done (%.1f sec elapsed)\n",
                  i, nrow(params_df), proc.time()["elapsed"] - t_batch_start))
    }
    sim_i <- safe_rxSolve_v3(model, get_p(row), ev, inits,
                              meta = list(drug = drug_name, N = n, seed = seed, id = row$id))
    if (is.null(sim_i)) return(NULL)
    sim_i$id <- row$id
    sim_i
  })
  sim_list <- Filter(Negate(is.null), sim_list)
  sim <- bind_rows(sim_list) %>% rename_with(~"conc", .cols = all_of(cfg$conc_col))

  c0_df <- sim %>% group_by(id) %>% slice(which.min(abs(time - t_disc))) %>%
    ungroup() %>% transmute(id, C0 = conc)
  post_df <- sim %>% filter(time >= t_disc) %>% mutate(post_time = time - t_disc) %>%
    select(id, post_time, conc) %>% left_join(c0_df, by = "id") %>%
    left_join(params_df %>% select(id, cyp2d6), by = "id")

  post_df %>%
    group_by(id, cyp2d6) %>%
    group_modify(~{
      times <- .x$post_time; concs <- .x$conc; C0 <- unique(.x$C0)
      r50 <- find_threshold_time_v3(times, concs, 0.50 * C0)
      r25 <- find_threshold_time_v3(times, concs, 0.25 * C0)
      tibble(t_50pct = r50$t, status_50pct = r50$status,
             status_25pct = r25$status)
    }) %>% ungroup()
}

# ── 3. N x seed 안정성 점검 ───────────────────────────────────────────────────
stab_rows <- list()
for (n in N_GRID) {
  for (seed in SEED_GRID) {
    for (drug in names(drug_config)) {
      cat(sprintf("[stability] %s | N=%d | seed=%d ...\n", drug, n, seed))
      t_combo_start <- proc.time()["elapsed"]
      res <- sample_and_run(drug, n, seed)
      cat(sprintf("  -> done (%.1f sec, %d patients skipped so far total)\n",
                  proc.time()["elapsed"] - t_combo_start, length(PATHOLOGICAL_LOG)))

      if (drug == "Paroxetine") {
        eff_r <- rank_biserial_r(
          res$t_50pct[res$cyp2d6 == "poor"   & res$status_50pct == "reached"],
          res$t_50pct[res$cyp2d6 == "normal" & res$status_50pct == "reached"]
        )
      } else {
        eff_r <- NA_real_
      }

      stab_rows[[length(stab_rows) + 1]] <- tibble(
        drug = drug, N = n, seed = seed,
        median_t50_h    = median(res$t_50pct[res$status_50pct == "reached"], na.rm = TRUE),
        pct_censored_25 = round(100 * mean(res$status_25pct == "censored", na.rm = TRUE), 2),
        pm_effect_size_r = eff_r
      )
    }
  }
}
stab_df <- bind_rows(stab_rows)

# 건너뛴(타임아웃/에러) 환자 기록 -- 있으면 반드시 파일로 남김(조용히 버리지 않음)
if (length(PATHOLOGICAL_LOG) > 0) {
  path_log_df <- bind_rows(lapply(PATHOLOGICAL_LOG, as_tibble))
  path_log_path <- file.path(OUT_DIR, "pathological_patients_v3.csv")
  write_csv(path_log_df, path_log_path)
  cat(sprintf("\n[WARN] 총 %d명의 환자가 %d초 제한시간을 넘겨 건너뛰어짐 -- %s에 기록함. 반환 보고서에 언급할 것.\n",
              length(PATHOLOGICAL_LOG), PATIENT_TIMEOUT_SEC, path_log_path))
} else {
  cat("\n모든 환자가 제한시간 안에 정상적으로 풀림 (건너뛴 환자 없음).\n")
}

cv_table <- stab_df %>%
  group_by(drug, N) %>%
  summarise(
    mean_median_t50  = round(mean(median_t50_h, na.rm = TRUE), 2),
    CV_t50_pct       = round(sd(median_t50_h, na.rm = TRUE) / mean(median_t50_h, na.rm = TRUE) * 100, 2),
    mean_pct_censored25 = round(mean(pct_censored_25, na.rm = TRUE), 2),
    sd_pct_censored25   = round(sd(pct_censored_25, na.rm = TRUE), 2),
    mean_pm_effect_r = round(mean(pm_effect_size_r, na.rm = TRUE), 3),
    sd_pm_effect_r   = round(sd(pm_effect_size_r, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(drug, N)

cat("\n=== Stability across N and seed (v3) ===\n")
print(cv_table, n = Inf, width = Inf)

raw_path <- file.path(OUT_DIR, "stability_raw_v3.csv")
cv_path  <- file.path(OUT_DIR, "stability_cv_v3.csv")
write_csv(stab_df,  raw_path)
write_csv(cv_table, cv_path)
cat(sprintf("\nSaved: %s\nSaved: %s\n", raw_path, cv_path))

# ── 4. Fluoxetine 보간 오차: dt=2h(주 분석 격자) vs dt=0.5h(세밀 격자) ───────
cat("\n=== Fluoxetine interpolation error: dt=2h vs dt=0.5h (sample n=20) ===\n")

set.seed(2025)
n_check <- 20
p <- pk_params$fluoxetine; nfx <- p$norfluoxetine
check_df <- tibble(
  id = seq_len(n_check),
  CL_sampled     = sample_lognormal_v3(n_check, p$CL, p$CV_CL),
  Vd_sampled     = sample_lognormal_v3(n_check, p$Vd, p$CV_Vd),
  CL_met_sampled = sample_lognormal_v3(n_check, nfx$CL_met, nfx$CV_CL_met),
  Vd_met_sampled = nfx$Vd_met
)

cfg <- drug_config[["Fluoxetine"]]
t_disc <- cfg$n_doses * cfg$tau_h

run_fluoxetine_dt <- function(dt_value) {
  ev <- et(amt = cfg$F * cfg$dose_mg, ii = cfg$tau_h, addl = cfg$n_doses - 1, time = 0) %>%
        et(seq(0, t_disc + cfg$obs_h, by = dt_value))
  sim_list <- lapply(seq_len(nrow(check_df)), function(i) {
    row <- check_df[i, , drop = FALSE]
    sim_i <- rxSolve(ode_parent_met,
                      c(Ka = cfg$Ka, Vd_p = row$Vd_sampled, CL_p = row$CL_sampled,
                        Fm = cfg$Fm, Vd_m = row$Vd_met_sampled, CL_m = row$CL_met_sampled),
                      ev, inits = c(depot = 0, Cp = 0, Cm = 0)) %>% as.data.frame()
    sim_i$id <- row$id
    sim_i
  })
  sim <- bind_rows(sim_list) %>% mutate(conc = Cp + Cm)

  c0_df <- sim %>% group_by(id) %>% slice(which.min(abs(time - t_disc))) %>%
    ungroup() %>% transmute(id, C0 = conc)
  post_df <- sim %>% filter(time >= t_disc) %>% mutate(post_time = time - t_disc) %>%
    select(id, post_time, conc) %>% left_join(c0_df, by = "id")

  post_df %>%
    group_by(id) %>%
    group_modify(~{
      times <- .x$post_time; concs <- .x$conc; C0 <- unique(.x$C0)
      r50 <- find_threshold_time_v3(times, concs, 0.50 * C0)
      tibble(C0 = C0, t_50pct = r50$t, status_50pct = r50$status)
    }) %>% ungroup()
}

res_dt2   <- run_fluoxetine_dt(2.0)  %>% rename(C0_dt2 = C0, t50_dt2 = t_50pct, status_dt2 = status_50pct)
res_dt05  <- run_fluoxetine_dt(0.5)  %>% rename(C0_dt05 = C0, t50_dt05 = t_50pct, status_dt05 = status_50pct)

interp_compare <- res_dt2 %>%
  left_join(res_dt05, by = "id") %>%
  mutate(
    C0_rel_diff_pct  = round(100 * (C0_dt2 - C0_dt05) / C0_dt05, 4),
    t50_abs_diff_h   = round(t50_dt2 - t50_dt05, 3),
    t50_rel_diff_pct = round(100 * (t50_dt2 - t50_dt05) / t50_dt05, 3)
  )

print(interp_compare, n = Inf, width = Inf)

interp_path <- file.path(OUT_DIR, "fluoxetine_interp_error_v3.csv")
write_csv(interp_compare, interp_path)
cat(sprintf("\nSaved: %s\n", interp_path))
cat(sprintf("Max |t50 rel diff|: %.3f%%  |  Max |C0 rel diff|: %.4f%%\n",
            max(abs(interp_compare$t50_rel_diff_pct), na.rm = TRUE),
            max(abs(interp_compare$C0_rel_diff_pct), na.rm = TRUE)))

cat("\nstability_v3.R 완료. 다음 단계: tables_v3.R\n\n")
