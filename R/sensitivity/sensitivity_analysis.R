# ============================================================
# Sensitivity analysis — fm,CYP2D6 가중치 + Potency-weighted active moiety
# 재시뮬레이션_필요항목_가이드.md 항목 3, 5 대응
# Author: Nari Kim (코드 작성 보조: Claude)
# ============================================================
#
# 목적:
#   항목 3, 5 모두 "정확한 문헌값을 못 찾음" → 범위 민감도분석으로 대체
#   (가이드 문서 89행: "가중치를 0.5~1.5 범위로 바꿔가며 결과가 얼마나
#    민감한지 보여주는 범위 민감도분석으로 대체해도 됨")
#
#   항목 3: 지금은 모든 약물에 동일한 CYP2D6 배수(PM×0.30/NM×1.00/RM×1.60)를
#           전체 CL에 곱함. fm,CYP2D6(그 약물 청소율 중 CYP2D6이 담당하는 비중)이
#           1.0(현재 가정)보다 작으면 결과가 얼마나 달라지는지 봄.
#   항목 5: 지금은 metabolite가 모약물과 동일한 역가라고 가정(potency=1.0).
#           0.5~1.5 범위에서 바뀌면 t50%C0 / 7일 AUC deficit이 얼마나 민감한지 봄.
#
# 사용법:
#   1. 아래 N_sens 값을 먼저 500 정도로 낮춰서 한 번 전체 실행해보고
#      (Source 버튼, 전체 스크립트 순서대로 실행 — 일부만 실행하지 말 것)
#      에러 없이 끝까지 돌면 N_sens를 2000으로 올려서 본 실행.
#   2. 끝나면 outputs/mc_sensitivity/ 아래 4개 CSV 생성됨.
#   3. item3_fm_cyp2d6_summary.csv, item5_potency_summary.csv 두 개를
#      (혹은 콘솔에 출력된 표 스크린샷을) 나한테 보내주면 원고에 반영할게.
#
# ============================================================

suppressPackageStartupMessages({
  library(rxode2)
  library(tidyverse)
})

BASE <- "~/nari-research/pkpd-antidepressant-sim"
dir.create(file.path(BASE, "outputs/mc_sensitivity"), showWarnings = FALSE, recursive = TRUE)

# 기존 repo의 PK 파라미터 재사용 (pk_params, cyp2d6 리스트)
source(file.path(BASE, "data/pk_parameters/pk_parameters.R"))

set.seed(2025)

# ------------------------------------------------------------
# 1. rxode2 모델 — 기존 mc_run_single.R의 ode_1cmt / ode_parent_met와
#    동일한 구조 (depot -> parent [-> metabolite], Ctotal = Cp+Cm)
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
# 2. drug_config — mc_run_single.R의 drug_config와 동일한 값
#    (paroxetine은 항목5 대상 아님 — 주요 활성대사체 없음)
# ------------------------------------------------------------
drug_config <- list(
  fluoxetine = list(
    model = "parent_met", Ka = 0.72, F = 0.70, Fm = 0.72,
    dose_mg = 20, tau_h = 24, n_doses = 60, obs_h = 60 * 24, dt = 1.0
  ),
  paroxetine = list(
    model = "1cmt", Ka = 0.580, F = 0.50,
    dose_mg = 20, tau_h = 24, n_doses = 14, obs_h = 14 * 24, dt = 0.5
  ),
  venlafaxine_IR = list(
    model = "parent_met", Ka = 1.50, F = 0.45, Fm = 0.80,
    dose_mg = 25, tau_h = 8, n_doses = 21, obs_h = 7 * 24, dt = 0.25
  ),
  venlafaxine_XR = list(
    model = "parent_met", Ka = 0.25, F = 0.45, Fm = 0.80,
    dose_mg = 75, tau_h = 24, n_doses = 7, obs_h = 7 * 24, dt = 0.25
  )
)

# ------------------------------------------------------------
# 3. 샘플링 함수 — mc_params.R과 동일한 로직
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

# 한 번만 샘플링 (CL_base = CYP2D6 배수 적용 전 clearance) -> 이후 fm_cyp2d6
# 값을 바꿔가며 같은 코호트에 다른 가중치만 적용 (불필요한 몬테카를로 노이즈 방지)
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

# 항목 3의 핵심 수정 지점:
# 기존: CL_sampled = CL_base * fac   (CYP2D6이 전체 CL의 100%를 좌우한다고 가정)
# 수정: CYP2D6이 담당하는 부분(fm_cyp2d6)만 배수 적용, 나머지(1-fm_cyp2d6)는 고정
apply_fm_cyp2d6 <- function(base_df, fm_cyp2d6) {
  base_df %>%
    mutate(CL_sampled = CL_base * ((1 - fm_cyp2d6) * 1.0 + fm_cyp2d6 * fac))
}

# ------------------------------------------------------------
# 4. 시뮬레이션 + 지표 추출 — mc_run_single.R의 run_batch와 동일한 로직
#    (find_threshold_time, trapz 포함). potency_weight만 새로 추가됨.
# ------------------------------------------------------------
find_threshold_time <- function(times, concs, threshold) {
  idx <- which(concs <= threshold)[1]
  if (is.na(idx) || idx == 1) return(NA_real_)
  t1 <- times[idx - 1]; t2 <- times[idx]
  c1 <- concs[idx - 1]; c2 <- concs[idx]
  t1 + (threshold - c1) * (t2 - t1) / (c2 - c1)
}

trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

# 항목 5의 핵심 수정 지점:
# 기존: conc = Cp + Cm              (metabolite가 모약물과 동등한 역가라고 가정)
# 수정: conc = Cp + potency_weight * Cm   (역가 가중치를 곱해서 sweep)
run_sensitivity_batch <- function(params_df, cfg, potency_weight = 1.0) {
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

  css_df <- sim %>%
    filter(time >= (cfg$n_doses - 1) * cfg$tau_h, time <= t_disc) %>%
    group_by(id) %>% summarise(Css = mean(conc, na.rm = TRUE), .groups = "drop")

  post_df <- sim %>%
    filter(time >= t_disc) %>%
    mutate(post_time = time - t_disc) %>%
    select(id, post_time, conc) %>%
    left_join(css_df, by = "id") %>%
    left_join(params_df %>% select(id, cyp2d6), by = "id")

  post_df %>%
    group_by(id, cyp2d6) %>%
    group_modify(~{
      times <- .x$post_time; concs <- .x$conc; Css <- unique(.x$Css)
      t50 <- find_threshold_time(times, concs, 0.5 * Css)
      post_7d <- .x %>% filter(post_time <= 7 * 24)
      auc_def <- if (nrow(post_7d) > 1) {
        Css * max(post_7d$post_time) - trapz(post_7d$post_time, post_7d$conc)
      } else {
        NA_real_
      }
      tibble(Css = Css, t50_C0 = t50, AUC_deficit_7d = auc_def)
    }) %>%
    ungroup()
}

# ============================================================
# 5. 항목 3 — fm,CYP2D6 sweep
# ============================================================
N_sens <- 500   # <- 처음엔 500으로 테스트, 에러 없으면 2000으로 올려서 재실행

fm_cyp2d6_grid <- c(0.3, 0.5, 0.7, 0.9, 1.0)  # 1.0 = 현재 원고의 가정(CL 100% CYP2D6-dependent)
drugs_cyp <- c("fluoxetine", "paroxetine", "venlafaxine_IR", "venlafaxine_XR")

item3_results <- list()
for (drug in drugs_cyp) {
  cfg <- drug_config[[drug]]
  base_df <- sample_drug_base(N_sens, drug)   # 같은 코호트를 fm 값마다 재사용
  for (fm in fm_cyp2d6_grid) {
    cat(sprintf("[항목3] %s | fm_cyp2d6=%.1f 시뮬레이션 중...\n", drug, fm))
    params_df <- apply_fm_cyp2d6(base_df, fm)
    res <- run_sensitivity_batch(params_df, cfg, potency_weight = 1.0)
    res$drug <- drug
    res$fm_cyp2d6 <- fm
    item3_results[[paste(drug, fm)]] <- res
  }
}
item3_all <- bind_rows(item3_results)

item3_summary <- item3_all %>%
  group_by(drug, fm_cyp2d6, cyp2d6) %>%
  summarise(median_t50_h = median(t50_C0, na.rm = TRUE), n = n(), .groups = "drop") %>%
  arrange(drug, cyp2d6, fm_cyp2d6)

cat("\n========== 항목 3: fm,CYP2D6 sweep 결과 (median t50%C0, h) ==========\n")
print(item3_summary, n = 100)

write.csv(item3_all,     file.path(BASE, "outputs/mc_sensitivity/item3_fm_cyp2d6_raw.csv"), row.names = FALSE)
write.csv(item3_summary, file.path(BASE, "outputs/mc_sensitivity/item3_fm_cyp2d6_summary.csv"), row.names = FALSE)

# ============================================================
# 6. 항목 5 — potency weight sweep
# ============================================================
potency_grid <- c(0.5, 0.75, 1.0, 1.25, 1.5)  # 1.0 = 현재 원고의 가정(동등 역가)
drugs_potency <- c("fluoxetine", "venlafaxine_IR", "venlafaxine_XR")

item5_results <- list()
for (drug in drugs_potency) {
  cfg <- drug_config[[drug]]
  base_df <- sample_drug_base(N_sens, drug)
  params_df <- apply_fm_cyp2d6(base_df, fm_cyp2d6 = 1.0)  # 현재 원고의 CYP2D6 가정 그대로 사용
  for (w in potency_grid) {
    cat(sprintf("[항목5] %s | potency_weight=%.2f 시뮬레이션 중...\n", drug, w))
    res <- run_sensitivity_batch(params_df, cfg, potency_weight = w)
    res$drug <- drug
    res$potency_weight <- w
    item5_results[[paste(drug, w)]] <- res
  }
}
item5_all <- bind_rows(item5_results)

item5_summary <- item5_all %>%
  group_by(drug, potency_weight) %>%
  summarise(
    median_t50_h  = median(t50_C0, na.rm = TRUE),
    median_AUCdef = median(AUC_deficit_7d, na.rm = TRUE),
    n = n(), .groups = "drop"
  ) %>%
  arrange(drug, potency_weight)

cat("\n========== 항목 5: potency weight sweep 결과 ==========\n")
print(item5_summary, n = 100)

write.csv(item5_all,     file.path(BASE, "outputs/mc_sensitivity/item5_potency_raw.csv"), row.names = FALSE)
write.csv(item5_summary, file.path(BASE, "outputs/mc_sensitivity/item5_potency_summary.csv"), row.names = FALSE)

cat("\n✅ 완료. outputs/mc_sensitivity/ 아래 4개 CSV 생성됨:\n")
cat("   item3_fm_cyp2d6_raw.csv / item3_fm_cyp2d6_summary.csv\n")
cat("   item5_potency_raw.csv / item5_potency_summary.csv\n")
cat("item3_fm_cyp2d6_summary.csv, item5_potency_summary.csv 두 개 (또는 콘솔 표 스크린샷) 보내주면 원고에 반영할게!\n")
