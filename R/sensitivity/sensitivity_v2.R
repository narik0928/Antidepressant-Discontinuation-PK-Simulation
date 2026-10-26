# ============================================================
# Sensitivity analysis v2 — fm,CYP2D6 가중치 + Potency-weighted active moiety
# Author: Nari Kim (코드 작성 보조: Claude)
#
# 원본 R/sensitivity/sensitivity_analysis.R은 건드리지 않음. 새 파일로 분리.
#
# v1 대비 바뀐 것 (HANDOFF_TO_NARI.md 6-3, 8-1 / GPT_STAGE2_COMPARISON.md A2)
#   A2 (핵심 수정) : v1의 AUC_deficit_7d는 절대값(mg/L·h)인데 원고 3.9절에는
#                    "normalised AUC deficit, %"로 잘못 표기돼 있었음. v2는
#                    절대값(AUC_deficit_7d_abs)과 환자별 정규화값
#                    (AUC_deficit_7d_norm_pct = 100×deficit/(C0×7일))을
#                    둘 다 분리해서 저장 — 더 이상 단위를 섞지 않음.
#   D1            : Css(마지막 투약구간 평균) → C0(중단 시점 순간 농도)로
#                    기준을 바꿈. mc_run_v2.R과 동일한 정의.
#   D3            : F를 drug_config에 하드코딩하지 않고 pk_parameters_v2.R의
#                    pk_params$<drug>$F를 그대로 읽음 (v1의 "F를 스크립트마다
#                    따로 적어서 생긴 하드코딩 문제" 재발 방지, 8-1 요구사항).
#                    v2에서는 F=1.0이므로 전체 공칭 용량이 들어감.
#   상태 코드      : find_threshold_time이 이제 reached/censored/calc_failed
#                    세 가지를 구분해서 반환 (v1은 둘 다 NA로 뭉뚱그림).
#   N 적정성 점검  : 환자 수(N_sens)가 충분한지 median에 대한 부트스트랩
#                    95% CI로 직접 확인 (요청: "N 적정성 판단용 bootstrap CI").
#
# 사용법: 1) N_sens를 500으로 먼저 전체 실행 → 에러 없이 끝까지 도는지 확인
#         2) item3_fm_cyp2d6_summary_v2.csv / item5_potency_summary_v2.csv의
#            CI 폭(ci_lo~ci_hi)을 보고 너무 넓으면 N_sens를 올려서 재실행
#            (가이드: 본 MC 코호트와 맞추려면 2000 정도)
# ============================================================

suppressPackageStartupMessages({
  library(rxode2)
  library(tidyverse)
})

BASE <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR <- file.path(BASE, "outputs", "mc_sensitivity_v2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "data/pk_parameters/pk_parameters_v2.R"))   # v2: F=1 통일본

set.seed(2025)

N_sens <- 500   # 먼저 500으로 테스트, 부트스트랩 CI 폭 보고 필요하면 2000으로

# ------------------------------------------------------------
# 1. rxode2 모델 (v1과 동일 구조)
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
# 2. drug_config — v2: F/Ka/Fm을 pk_parameters_v2.R에서 직접 읽음
#    (하드코딩 금지, 8-1 요구사항). dose_mg/tau_h/n_doses는 숫자 유지(D2).
#    Fluoxetine만 D4(i)에 따라 obs_h/dt를 mc_run_v2.R과 동일하게 확장.
# ------------------------------------------------------------
drug_config <- list(
  fluoxetine = list(
    model = "parent_met",
    Ka = pk_params$fluoxetine$Ka,
    F  = pk_params$fluoxetine$F,                     # v2: 1.0
    Fm = pk_params$fluoxetine$norfluoxetine$Fm,
    dose_mg = 20, tau_h = 24, n_doses = 60,
    obs_h = 400 * 24, dt = 2.0                        # v2 (D4-i): mc_run_v2.R과 통일
  ),
  paroxetine = list(
    model = "1cmt",
    Ka = pk_params$paroxetine$Ka,
    F  = pk_params$paroxetine$F,                      # v2: 1.0
    dose_mg = 20, tau_h = 24, n_doses = 14, obs_h = 14 * 24, dt = 0.5
  ),
  venlafaxine_IR = list(
    model = "parent_met",
    Ka = pk_params$venlafaxine_IR$Ka,
    F  = pk_params$venlafaxine_IR$F,                  # v2: 1.0
    Fm = pk_params$venlafaxine_IR$ODV$Fm,
    dose_mg = 25, tau_h = 8, n_doses = 21, obs_h = 7 * 24, dt = 0.25
  ),
  venlafaxine_XR = list(
    model = "parent_met",
    Ka = pk_params$venlafaxine_XR$Ka,
    F  = pk_params$venlafaxine_XR$F,                  # v2: 1.0
    Fm = pk_params$venlafaxine_XR$ODV$Fm,
    dose_mg = 75, tau_h = 24, n_doses = 7, obs_h = 7 * 24, dt = 0.25
  )
)

# ------------------------------------------------------------
# 3. 샘플링 함수 (v1과 동일 로직 — mc_params.R 기반)
# ------------------------------------------------------------
assign_cyp2d6 <- function(n) {
  sample(x = c("poor", "normal", "rapid"), size = n, replace = TRUE,
         prob = c(cyp2d6$poor_metabolizer$freq,
                  cyp2d6$normal_metabolizer$freq,
                  cyp2d6$rapid_metabolizer$freq))
}

get_cyp2d6_factor <- function(phenotype) {
  case_when(
    phenotype == "poor"   ~ cyp2d6$poor_metabolizer$CL_factor,
    phenotype == "normal" ~ cyp2d6$normal_metabolizer$CL_factor,
    phenotype == "rapid"  ~ cyp2d6$rapid_metabolizer$CL_factor
  )
}

sample_lognormal <- function(n, pop_value, cv) {
  sigma <- sqrt(log(cv^2 + 1))
  mu <- log(pop_value) - sigma^2 / 2
  rlnorm(n, meanlog = mu, sdlog = sigma)
}

sample_drug_base <- function(n, drug_key) {
  p   <- pk_params[[drug_key]]
  cyp <- assign_cyp2d6(n)
  fac <- get_cyp2d6_factor(cyp)

  base <- tibble(
    id = seq_len(n), cyp2d6 = cyp, fac = fac,
    CL_base    = sample_lognormal(n, p$CL, p$CV_CL),
    Vd_sampled = sample_lognormal(n, p$Vd, p$CV_Vd)
  )

  if (drug_key == "fluoxetine") {
    met <- p$norfluoxetine
    base <- base %>% mutate(
      CL_met_sampled = sample_lognormal(n, met$CL_met, met$CV_CL_met),
      Vd_met_sampled = rep(met$Vd_met, n)
    )
  }
  if (drug_key %in% c("venlafaxine_IR", "venlafaxine_XR")) {
    met <- p$ODV
    base <- base %>% mutate(
      CL_met_sampled = sample_lognormal(n, met$CL_odv, met$CV_CL_odv),
      Vd_met_sampled = rep(met$Vd_odv, n)
    )
  }
  base
}

apply_fm_cyp2d6 <- function(base_df, fm_cyp2d6) {
  base_df %>%
    mutate(CL_sampled = CL_base * ((1 - fm_cyp2d6) * 1.0 + fm_cyp2d6 * fac))
}

# ------------------------------------------------------------
# 4. 시뮬레이션 + 지표 추출 (v2: C0 기준 + 상태 코드 + AUC 정규화 분리)
# ------------------------------------------------------------

# v2: idx==1(이론상 발생 불가)과 "끝까지 미도달"(censored)을 구분
find_threshold_time_v2 <- function(times, concs, threshold) {
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

# 항목 5 로직 동일: conc = Cp + potency_weight * Cm
run_sensitivity_batch_v2 <- function(params_df, cfg, potency_weight = 1.0) {
  t_disc <- cfg$n_doses * cfg$tau_h
  ev <- et(amt = cfg$F * cfg$dose_mg, ii = cfg$tau_h, addl = cfg$n_doses - 1, time = 0) %>%
        et(seq(0, t_disc + cfg$obs_h, by = cfg$dt))

  if (cfg$model == "1cmt") {
    model <- ode_1cmt
    inits <- c(depot = 0, A_c = 0)
    get_p <- function(row) c(Ka = cfg$Ka, Vd = row$Vd_sampled, CL = row$CL_sampled)
  } else {
    model <- ode_parent_met
    inits <- c(depot = 0, A_p = 0, A_m = 0)
    get_p <- function(row) c(Ka = cfg$Ka, Vd_p = row$Vd_sampled, CL_p = row$CL_sampled,
                              Fm = cfg$Fm, Vd_m = row$Vd_met_sampled, CL_m = row$CL_met_sampled)
  }

  sim_list <- lapply(seq_len(nrow(params_df)), function(i) {
    row <- params_df[i, , drop = FALSE]
    sim_i <- rxSolve(model, get_p(row), ev, inits = inits) %>% as.data.frame()
    sim_i$id <- row$id
    sim_i
  })
  sim <- bind_rows(sim_list)

  sim <- if (cfg$model == "1cmt") {
    sim %>% mutate(conc = Cp)
  } else {
    sim %>% mutate(conc = Cp + potency_weight * Cm)
  }

  # v2: C0 = 중단 시점(t_disc)과 가장 가까운 격자점의 농도 (mc_run_v2.R과 동일 정의)
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
    left_join(params_df %>% select(id, cyp2d6), by = "id")

  obs_7d_cap <- min(168, cfg$obs_h)

  post_df %>%
    group_by(id, cyp2d6) %>%
    group_modify(~{
      times <- .x$post_time; concs <- .x$conc; C0 <- unique(.x$C0)

      r50 <- find_threshold_time_v2(times, concs, 0.5 * C0)

      obs_7d  <- min(obs_7d_cap, max(times, na.rm = TRUE))
      post_7d <- .x %>% filter(post_time <= obs_7d)

      if (nrow(post_7d) >= 2 && !is.na(C0) && C0 > 0) {
        auc_actual   <- trapz(post_7d$post_time, post_7d$conc)
        auc_expected <- C0 * obs_7d
        auc_def_abs  <- auc_expected - auc_actual
        auc_def_norm <- 100 * auc_def_abs / auc_expected     # v2 (A2 fix): 환자별 정규화 %
      } else {
        auc_def_abs  <- NA_real_
        auc_def_norm <- NA_real_
      }

      tibble(
        C0 = C0,
        t50_C0 = r50$t, status_t50 = r50$status,
        AUC_deficit_7d_abs = auc_def_abs,            # mg/L·h — v1의 "AUC_deficit_7d"와 동일 단위(절대값)
        AUC_deficit_7d_norm_pct = auc_def_norm        # % — v1이 잘못 "정규화"라고 적었던 실제 값
      )
    }) %>%
    ungroup()
}

# 중앙값에 대한 부트스트랩 95% CI (N_sens 적정성 판단용)
bootstrap_median_ci <- function(x, B = 2000, conf = 0.95) {
  x <- x[!is.na(x)]
  if (length(x) < 5) return(c(lo = NA_real_, hi = NA_real_))
  boots <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  alpha <- (1 - conf) / 2
  c(lo = unname(quantile(boots, alpha)), hi = unname(quantile(boots, 1 - alpha)))
}

# ============================================================
# 5. 항목 3 — fm,CYP2D6 sweep (v2)
# ============================================================
fm_cyp2d6_grid <- c(0.3, 0.5, 0.7, 0.9, 1.0)
drugs_cyp <- c("fluoxetine", "paroxetine", "venlafaxine_IR", "venlafaxine_XR")

item3_results <- list()
for (drug in drugs_cyp) {
  cfg <- drug_config[[drug]]
  base_df <- sample_drug_base(N_sens, drug)
  for (fm in fm_cyp2d6_grid) {
    cat(sprintf("[v2 항목3] %s | fm_cyp2d6=%.1f 시뮬레이션 중...\n", drug, fm))
    params_df <- apply_fm_cyp2d6(base_df, fm)
    res <- run_sensitivity_batch_v2(params_df, cfg, potency_weight = 1.0)
    res$drug <- drug
    res$fm_cyp2d6 <- fm
    item3_results[[paste(drug, fm)]] <- res
  }
}
item3_all <- bind_rows(item3_results)

item3_summary <- item3_all %>%
  group_by(drug, fm_cyp2d6, cyp2d6) %>%
  summarise(
    n_total       = n(),
    n_reached     = sum(status_t50 == "reached", na.rm = TRUE),
    n_censored    = sum(status_t50 == "censored", na.rm = TRUE),
    n_calc_failed = sum(status_t50 == "calc_failed", na.rm = TRUE),
    median_t50_h  = median(t50_C0[status_t50 == "reached"], na.rm = TRUE),
    t50_ci_lo     = bootstrap_median_ci(t50_C0[status_t50 == "reached"])["lo"],
    t50_ci_hi     = bootstrap_median_ci(t50_C0[status_t50 == "reached"])["hi"],
    .groups = "drop"
  ) %>%
  arrange(drug, cyp2d6, fm_cyp2d6)

cat("\n========== v2 항목 3: fm,CYP2D6 sweep 결과 (median t50%C0, h + bootstrap 95% CI) ==========\n")
print(item3_summary, n = 100)

write.csv(item3_all,     file.path(OUT_DIR, "item3_fm_cyp2d6_raw_v2.csv"), row.names = FALSE)
write.csv(item3_summary, file.path(OUT_DIR, "item3_fm_cyp2d6_summary_v2.csv"), row.names = FALSE)

# ============================================================
# 6. 항목 5 — potency weight sweep (v2)
# ============================================================
potency_grid <- c(0.5, 0.75, 1.0, 1.25, 1.5)
drugs_potency <- c("fluoxetine", "venlafaxine_IR", "venlafaxine_XR")

item5_results <- list()
for (drug in drugs_potency) {
  cfg <- drug_config[[drug]]
  base_df <- sample_drug_base(N_sens, drug)
  params_df <- apply_fm_cyp2d6(base_df, fm_cyp2d6 = 1.0)
  for (w in potency_grid) {
    cat(sprintf("[v2 항목5] %s | potency_weight=%.2f 시뮬레이션 중...\n", drug, w))
    res <- run_sensitivity_batch_v2(params_df, cfg, potency_weight = w)
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
    median_AUCdef_norm_pct = median(AUC_deficit_7d_norm_pct, na.rm = TRUE),   # v2: 정규화 %값 (A2 fix)
    AUCdef_norm_ci_lo = bootstrap_median_ci(AUC_deficit_7d_norm_pct)["lo"],
    AUCdef_norm_ci_hi = bootstrap_median_ci(AUC_deficit_7d_norm_pct)["hi"],
    median_AUCdef_abs = median(AUC_deficit_7d_abs, na.rm = TRUE),             # mg/L·h, 참고용
    .groups = "drop"
  ) %>%
  arrange(drug, potency_weight)

cat("\n========== v2 항목 5: potency weight sweep 결과 (AUC deficit은 %와 mg/L·h 분리 표기) ==========\n")
print(item5_summary, n = 100)

write.csv(item5_all,     file.path(OUT_DIR, "item5_potency_raw_v2.csv"), row.names = FALSE)
write.csv(item5_summary, file.path(OUT_DIR, "item5_potency_summary_v2.csv"), row.names = FALSE)

# ============================================================
# 7. N 적정성 메모
# ============================================================
cat("\n── N_sens 적정성 점검 ──\n")
cat(sprintf("현재 N_sens = %d. 위 두 요약표의 ci_lo~ci_hi 폭이 median 대비\n", N_sens))
cat("넓게 보이면(특히 censored 비율이 높은 phenotype/weight 조합) N_sens를\n")
cat("2000으로 올려서 다시 실행한 뒤 CI 폭이 좁아지는지 확인.\n")

cat("\n✅ 완료. outputs/mc_sensitivity_v2/ 아래 4개 CSV 생성됨.\n")
cat("item3_fm_cyp2d6_summary_v2.csv, item5_potency_summary_v2.csv 두 개 보내주면\n")
cat("원고 3.9절/Limitation을 v2 숫자로 다시 쓸게 (단, 전체 재실행 결과가 먼저 나온 뒤에).\n")
