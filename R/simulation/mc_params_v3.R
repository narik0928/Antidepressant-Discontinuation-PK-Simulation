# =============================================================================
# mc_params_v3.R
# Phase 4 (v3) -- Monte Carlo Parameter Re-sampling with Z-value Recovery
# Author: Nari Kim (code assist: Claude)
#
# HANDOFF_TO_NARI_v3.md 7-1-2 요구사항 구현:
#   "같은 가상 환자를 유지합니다. 파라미터 중앙값과 CV가 바뀌었으므로,
#   v1.v2와 같은 표준정규 난수(z값)를 새 뮤, 시그마에 다시 적용하세요.
#   이 방식을 Methods에 한 문장으로 적습니다. 기존 파일에 z값이 없어
#   복원할 수 없으면 seed 2025로 새로 뽑고 그렇게 적습니다."
#   "로그정규 샘플링을 중앙값보존 뮤=log(population), 시그마^2=log(1+CV^2)로
#   바꿉니다(D11). v2의 평균보존식(-시그마^2/2)은 쓰지 않습니다."
#   "CYP2D6: paroxetine만 PM x 0.58. fluoxetine과 venlafaxine에는 phenotype을
#   배정하지 않습니다."
#
# 구조 메모: v2까지는 mc_params.R(샘플링, 1회 실행 후 결과를 mc_params_all_
# n10000_v1.csv로 동결)과 mc_run_v2.R(그 동결 파일을 읽어 ODE만 실행)이
# 분리돼 있었다. v3는 모집단값.CV.CYP2D6 구조 자체가 바뀌므로 반드시 다시
# 샘플링해야 하지만, "같은 가상 환자 유지(z값 재적용)"는 순수한 데이터
# 변환이라 ODE 실행(mc_run_v3.R)과 분리해 이 파일에서 처리한다.
# mc_run_v3.R은 이 파일을 source()하거나 이 파일의 출력 CSV를 읽는다
# (mc_run_v2.R이 mc_params_all_n10000_v1.csv를 읽던 것과 같은 패턴).
#
# z값 복원 원리:
#   mc_params.R의 sample_lognormal()은 rlnorm(n, meanlog=mu, sdlog=sigma)를
#   썼다. 이는 exp(mu + sigma*z) (z ~ N(0,1))와 수학적으로 같다. v1 산출물
#   (outputs/mc/mc_params_all_n10000_v1.csv)에는 z 자체는 저장돼 있지 않지만,
#   CL_sampled/Vd_sampled(이미 뽑힌 값)와 CL_pop/Vd_pop(당시 모집단값)이
#   전부 남아 있고, 당시 CV와 CYP2D6 배수는 pk_parameters_v2.R에 그대로
#   남아 있다(D2: "CL.Vd 등 숫자 값 자체는 전부 유지" -- v1과 v2는 수치가
#   동일함, 아래 STARTUP CHECK에서 재확인). 그래서
#     z = (log(raw_draw) - mu_old) / sigma_old
#   로 역산할 수 있다(raw_draw는 CYP2D6 배수가 곱해져 있었다면 먼저 그
#   배수로 나눠 되돌린 값). 이렇게 복원한 z를 새 뮤'(D11: 중앙값보존,
#   = log(pop_v3)), 시그마'(= sqrt(log(1+CV_v3^2)))에 그대로 적용하면
#   "같은 환자"가 유지된다(exp(mu_old+sigma_old*z) = raw_draw 이므로
#   z 자체는 CV·모집단값이 바뀌어도 "그 환자가 분포에서 몇 표준편차
#   위치에 있었는지"를 그대로 보존한다).
#
# CYP2D6 phenotype 재배정 (이 파일의 설계 판단 -- WORK_LOG_v3.md에도 기록.
# HANDOFF_v3가 직접 지정한 절차는 아니나 "같은 가상 환자 유지" 원칙과
# 가장 합리적으로 들어맞는 방법이라고 판단함. 7-2 전형환자 대조 후
# 덕현과 공유할 것):
#   - sertraline: 변경 없음 (원래도 CYP2D6 미적용, cyp2d6 = "normal" 고정)
#   - paroxetine: v1/v2에서 이미 뽑힌 phenotype 중 "poor"는 그대로
#     "poor"(PM)로 유지하고, "normal"과 "rapid"는 모두 "normal"(non-PM)로
#     합친다. PM 빈도(8%)가 v2와 v3에서 숫자가 같으므로, "이 환자가
#     PM이었는가"라는 이진 정체성이 그대로 보존된다. RM 범주 자체가
#     삭제됐으므로(D5) 그 구분만 없앤다.
#   - fluoxetine, venlafaxine_IR, venlafaxine_XR: D5에 따라 CYP2D6를 더
#     이상 모델링하지 않으므로, 옛 phenotype 라벨이 무엇이었든 배수는
#     전부 1.0으로 고정한다. 옛 라벨은 cyp2d6_v1_label 컬럼에 참고용으로만
#     남긴다(예: paroxetine을 제외한 나머지 약물에서 phenotype 효과가
#     실제로 사라졌는지 재현성 점검에 쓸 수 있음).
#
# 복원 실패 시(outputs/mc/mc_params_all_n10000_v1.csv가 없으면): 아래
# FRESH SAMPLE 경로로 전환해 seed 2025로 새로 뽑고, 콘솔과 산출 CSV의
# sampling_mode 컬럼에 "fresh_seed2025"라고 명시한다(HANDOFF_v3 7-1-2).
#
# Output:
#   outputs/mc_v3/mc_params_v3_<drug>_n<N>.csv
#   outputs/mc_v3/mc_params_v3_all_n<N>.csv
#
# 실행: source("R/simulation/mc_params_v3.R")
# 먼저 N=100으로 에러 없이 끝까지 도는지 확인한 뒤 N=10000으로.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

BASE    <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR <- file.path(BASE, "outputs", "mc_v3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N <- 100   # 먼저 100으로 테스트 (7-2), 통과하면 10000으로 바꿔서 본 실행

RECOVERY_FILE <- file.path(BASE, "outputs", "mc", "mc_params_all_n10000_v1.csv")

# 새 파라미터(v3) -- 전역 pk_params, cyp2d6 로 들어옴
source(file.path(BASE, "data", "pk_parameters", "pk_parameters_v3.R"))
pk_params_v3 <- pk_params
cyp2d6_v3    <- cyp2d6
rm(pk_params, cyp2d6)

# 역산에 필요한 "옛" 파라미터(v2 -- D2 결정상 v1과 숫자 동일) -- 전역을
# 덮어쓰지 않도록 별도 환경에 로드한다.
old_env <- new.env()
sys.source(file.path(BASE, "data", "pk_parameters", "pk_parameters_v2.R"), envir = old_env)

# ── STARTUP CHECK: v1 CSV의 CL_pop/Vd_pop이 pk_parameters_v2.R의 값과
#    실제로 일치하는지 확인 (D2 "숫자 유지" 전제가 깨졌으면 복원 자체가
#    틀어지므로 반드시 먼저 확인) ─────────────────────────────────────────
verify_v1_consistency <- function(old_all) {
  checks <- list(
    list(drug = "Sertraline",      p = old_env$pk_params$sertraline),
    list(drug = "Fluoxetine",      p = old_env$pk_params$fluoxetine),
    list(drug = "Paroxetine",      p = old_env$pk_params$paroxetine),
    list(drug = "Venlafaxine_IR",  p = old_env$pk_params$venlafaxine_IR),
    list(drug = "Venlafaxine_XR",  p = old_env$pk_params$venlafaxine_XR)
  )
  ok <- TRUE
  for (chk in checks) {
    row1 <- old_all %>% filter(drug == chk$drug) %>% slice(1)
    if (nrow(row1) == 0) next
    if (abs(row1$CL_pop - chk$p$CL) > 1e-6 || abs(row1$Vd_pop - chk$p$Vd) > 1e-6) {
      warning(sprintf(
        "[mc_params_v3] %s: v1 CSV의 CL_pop/Vd_pop(%.4f/%.4f)이 pk_parameters_v2.R(%.4f/%.4f)과 다릅니다 -- z 복원이 맞지 않을 수 있습니다. 직접 확인하세요.",
        chk$drug, row1$CL_pop, row1$Vd_pop, chk$p$CL, chk$p$Vd
      ))
      ok <- FALSE
    }
  }
  if (ok) cat("STARTUP CHECK 통과: v1 CSV의 CL_pop/Vd_pop이 pk_parameters_v2.R과 일치합니다.\n")
  ok
}

# ── 로그정규 변환 헬퍼 ────────────────────────────────────────────────────────

# v1/v2 방식(평균보존): mu = log(pop) - sigma^2/2
old_mu_sigma <- function(pop_value, cv) {
  sigma <- sqrt(log(cv^2 + 1))
  mu    <- log(pop_value) - sigma^2 / 2
  list(mu = mu, sigma = sigma)
}

# v3 방식(중앙값보존, D11): mu = log(pop)
new_mu_sigma <- function(pop_value, cv) {
  sigma <- sqrt(log(cv^2 + 1))
  mu    <- log(pop_value)
  list(mu = mu, sigma = sigma)
}

recover_z <- function(sampled_value, pop_value_old, cv_old) {
  ms <- old_mu_sigma(pop_value_old, cv_old)
  (log(sampled_value) - ms$mu) / ms$sigma
}

reapply_z <- function(z, pop_value_new, cv_new) {
  ms <- new_mu_sigma(pop_value_new, cv_new)
  exp(ms$mu + ms$sigma * z)
}

sample_lognormal_v3 <- function(n, pop_value, cv) {
  ms <- new_mu_sigma(pop_value, cv)
  rlnorm(n, meanlog = ms$mu, sdlog = ms$sigma)
}

# 옛 CYP2D6 배수(v1/v2, 3범주: poor/normal/rapid) -- 역산용
old_cyp_factor <- function(phenotype) {
  case_when(
    phenotype == "poor"   ~ old_env$cyp2d6$poor_metabolizer$CL_factor,   # 0.30
    phenotype == "normal" ~ old_env$cyp2d6$normal_metabolizer$CL_factor, # 1.00
    phenotype == "rapid"  ~ old_env$cyp2d6$rapid_metabolizer$CL_factor,  # 1.60
    TRUE ~ 1.0
  )
}

# ── 2. 약물별 복원 함수 (RECOVERED 경로) ──────────────────────────────────────

recover_sertraline_v3 <- function(old_all) {
  df <- old_all %>% filter(drug == "Sertraline")
  p_old <- old_env$pk_params$sertraline
  p_new <- pk_params_v3$sertraline

  z_CL <- recover_z(df$CL_sampled, p_old$CL, p_old$CV_CL)
  z_Vd <- recover_z(df$Vd_sampled, p_old$Vd, p_old$CV_Vd)

  tibble(
    id              = df$id,
    drug            = "Sertraline",
    cyp2d6          = "normal",         # CYP2D6 미적용, v2와 동일하게 고정
    cyp2d6_v1_label = df$cyp2d6,
    CL_pop          = p_new$CL,
    Vd_pop          = p_new$Vd,
    CL_sampled      = reapply_z(z_CL, p_new$CL, p_new$CV_CL),
    Vd_sampled      = reapply_z(z_Vd, p_new$Vd, p_new$CV_Vd),
    CL_met_sampled  = NA_real_,
    Vd_met_sampled  = NA_real_,
    z_CL = z_CL, z_Vd = z_Vd, z_CL_met = NA_real_
  )
}

recover_fluoxetine_v3 <- function(old_all) {
  df <- old_all %>% filter(drug == "Fluoxetine")
  p_old   <- old_env$pk_params$fluoxetine
  nfx_old <- old_env$pk_params$fluoxetine$norfluoxetine
  p_new   <- pk_params_v3$fluoxetine
  nfx_new <- pk_params_v3$fluoxetine$norfluoxetine

  old_fac <- old_cyp_factor(df$cyp2d6)          # v1/v2: CYP2D6 -> parent CL에 적용
  raw_CL  <- df$CL_sampled / old_fac
  z_CL    <- recover_z(raw_CL, p_old$CL, p_old$CV_CL)
  z_Vd    <- recover_z(df$Vd_sampled, p_old$Vd, p_old$CV_Vd)
  z_CLmet <- recover_z(df$CL_met_sampled, nfx_old$CL_met, nfx_old$CV_CL_met)  # 배수 없음

  tibble(
    id              = df$id,
    drug            = "Fluoxetine",
    cyp2d6          = "not_modeled",    # D5: fluoxetine은 phenotype 배정 안 함
    cyp2d6_v1_label = df$cyp2d6,        # 참고용(과거 라벨)
    CL_pop          = p_new$CL,
    Vd_pop          = p_new$Vd,
    CL_sampled      = reapply_z(z_CL, p_new$CL, p_new$CV_CL),      # x 1.0 (D5, 배수 없음)
    Vd_sampled      = reapply_z(z_Vd, p_new$Vd, p_new$CV_Vd),
    CL_met_sampled  = reapply_z(z_CLmet, nfx_new$CL_met, nfx_new$CV_CL_met),
    Vd_met_sampled  = nfx_new$Vd_met,    # 고정값(IIV 없음), v3 값(값 유지, 2310)
    z_CL = z_CL, z_Vd = z_Vd, z_CL_met = z_CLmet
  )
}

recover_paroxetine_v3 <- function(old_all) {
  df <- old_all %>% filter(drug == "Paroxetine")
  p_old <- old_env$pk_params$paroxetine
  p_new <- pk_params_v3$paroxetine

  old_fac <- old_cyp_factor(df$cyp2d6)
  raw_CL  <- df$CL_sampled / old_fac
  z_CL    <- recover_z(raw_CL, p_old$CL, p_old$CV_CL)
  z_Vd    <- recover_z(df$Vd_sampled, p_old$Vd, p_old$CV_Vd)

  # phenotype 재배정: poor -> poor(PM) 유지, normal/rapid -> normal(non-PM) 합침
  cyp_v3  <- if_else(df$cyp2d6 == "poor", "poor", "normal")
  new_fac <- if_else(cyp_v3 == "poor",
                      cyp2d6_v3$poor_metabolizer$CL_factor,     # 0.58 (D5, Sindrup 1992)
                      cyp2d6_v3$normal_metabolizer$CL_factor)   # 1.00

  tibble(
    id              = df$id,
    drug            = "Paroxetine",
    cyp2d6          = cyp_v3,
    cyp2d6_v1_label = df$cyp2d6,
    CL_pop          = p_new$CL,
    Vd_pop          = p_new$Vd,
    CL_sampled      = reapply_z(z_CL, p_new$CL, p_new$CV_CL) * new_fac,
    Vd_sampled      = reapply_z(z_Vd, p_new$Vd, p_new$CV_Vd),
    CL_met_sampled  = NA_real_,
    Vd_met_sampled  = NA_real_,
    z_CL = z_CL, z_Vd = z_Vd, z_CL_met = NA_real_
  )
}

recover_venlafaxine_v3 <- function(old_all, formulation) {
  # formulation: "IR" 또는 "XR"
  drug_name <- paste0("Venlafaxine_", formulation)
  df <- old_all %>% filter(drug == drug_name)
  p_old <- if (formulation == "IR") old_env$pk_params$venlafaxine_IR else old_env$pk_params$venlafaxine_XR
  p_new <- if (formulation == "IR") pk_params_v3$venlafaxine_IR      else pk_params_v3$venlafaxine_XR
  odv_old <- p_old$ODV
  odv_new <- p_new$ODV

  old_fac  <- old_cyp_factor(df$cyp2d6)      # v1/v2: CYP2D6 -> parent CL에 적용
  raw_CL   <- df$CL_sampled / old_fac
  z_CL     <- recover_z(raw_CL, p_old$CL, p_old$CV_CL)
  z_Vd     <- recover_z(df$Vd_sampled, p_old$Vd, p_old$CV_Vd)
  z_CLodv  <- recover_z(df$CL_met_sampled, odv_old$CL_odv, odv_old$CV_CL_odv)  # 배수 없음

  tibble(
    id              = df$id,
    drug            = drug_name,
    cyp2d6          = "not_modeled",    # D5: venlafaxine은 phenotype 배정 안 함
    cyp2d6_v1_label = df$cyp2d6,
    CL_pop          = p_new$CL,
    Vd_pop          = p_new$Vd,
    CL_sampled      = reapply_z(z_CL, p_new$CL, p_new$CV_CL),      # x 1.0 (D5)
    Vd_sampled      = reapply_z(z_Vd, p_new$Vd, p_new$CV_Vd),
    CL_met_sampled  = reapply_z(z_CLodv, odv_new$CL_odv, odv_new$CV_CL_odv),
    Vd_met_sampled  = odv_new$Vd_odv,    # 고정값(IIV 없음), v3 값(값 유지, 399)
    z_CL = z_CL, z_Vd = z_Vd, z_CL_met = z_CLodv
  )
}

# ── 3. 신규 샘플링 함수 (FRESH SAMPLE 경로, 복원 불가 시) ────────────────────

assign_cyp2d6_v3_paroxetine <- function(n) {
  sample(
    x       = c("poor", "normal"),
    size    = n,
    replace = TRUE,
    prob    = c(cyp2d6_v3$poor_metabolizer$freq, cyp2d6_v3$normal_metabolizer$freq)
  )
}

fresh_sample_all_v3 <- function(n) {
  cat("FRESH SAMPLE (recovery file not found) -- seed 2025로 새로 뽑습니다.\n")
  set.seed(2025)

  p <- pk_params_v3

  sertraline <- tibble(
    id = seq_len(n), drug = "Sertraline", cyp2d6 = "normal", cyp2d6_v1_label = NA_character_,
    CL_pop = p$sertraline$CL, Vd_pop = p$sertraline$Vd,
    CL_sampled = sample_lognormal_v3(n, p$sertraline$CL, p$sertraline$CV_CL),
    Vd_sampled = sample_lognormal_v3(n, p$sertraline$Vd, p$sertraline$CV_Vd),
    CL_met_sampled = NA_real_, Vd_met_sampled = NA_real_,
    z_CL = NA_real_, z_Vd = NA_real_, z_CL_met = NA_real_
  )

  fluoxetine <- tibble(
    id = seq_len(n), drug = "Fluoxetine", cyp2d6 = "not_modeled", cyp2d6_v1_label = NA_character_,
    CL_pop = p$fluoxetine$CL, Vd_pop = p$fluoxetine$Vd,
    CL_sampled = sample_lognormal_v3(n, p$fluoxetine$CL, p$fluoxetine$CV_CL),
    Vd_sampled = sample_lognormal_v3(n, p$fluoxetine$Vd, p$fluoxetine$CV_Vd),
    CL_met_sampled = sample_lognormal_v3(n, p$fluoxetine$norfluoxetine$CL_met, p$fluoxetine$norfluoxetine$CV_CL_met),
    Vd_met_sampled = p$fluoxetine$norfluoxetine$Vd_met,
    z_CL = NA_real_, z_Vd = NA_real_, z_CL_met = NA_real_
  )

  cyp_par <- assign_cyp2d6_v3_paroxetine(n)
  fac_par <- if_else(cyp_par == "poor", cyp2d6_v3$poor_metabolizer$CL_factor, cyp2d6_v3$normal_metabolizer$CL_factor)
  paroxetine <- tibble(
    id = seq_len(n), drug = "Paroxetine", cyp2d6 = cyp_par, cyp2d6_v1_label = NA_character_,
    CL_pop = p$paroxetine$CL, Vd_pop = p$paroxetine$Vd,
    CL_sampled = sample_lognormal_v3(n, p$paroxetine$CL, p$paroxetine$CV_CL) * fac_par,
    Vd_sampled = sample_lognormal_v3(n, p$paroxetine$Vd, p$paroxetine$CV_Vd),
    CL_met_sampled = NA_real_, Vd_met_sampled = NA_real_,
    z_CL = NA_real_, z_Vd = NA_real_, z_CL_met = NA_real_
  )

  ven_ir <- tibble(
    id = seq_len(n), drug = "Venlafaxine_IR", cyp2d6 = "not_modeled", cyp2d6_v1_label = NA_character_,
    CL_pop = p$venlafaxine_IR$CL, Vd_pop = p$venlafaxine_IR$Vd,
    CL_sampled = sample_lognormal_v3(n, p$venlafaxine_IR$CL, p$venlafaxine_IR$CV_CL),
    Vd_sampled = sample_lognormal_v3(n, p$venlafaxine_IR$Vd, p$venlafaxine_IR$CV_Vd),
    CL_met_sampled = sample_lognormal_v3(n, p$venlafaxine_IR$ODV$CL_odv, p$venlafaxine_IR$ODV$CV_CL_odv),
    Vd_met_sampled = p$venlafaxine_IR$ODV$Vd_odv,
    z_CL = NA_real_, z_Vd = NA_real_, z_CL_met = NA_real_
  )

  ven_xr <- tibble(
    id = seq_len(n), drug = "Venlafaxine_XR", cyp2d6 = "not_modeled", cyp2d6_v1_label = NA_character_,
    CL_pop = p$venlafaxine_XR$CL, Vd_pop = p$venlafaxine_XR$Vd,
    CL_sampled = sample_lognormal_v3(n, p$venlafaxine_XR$CL, p$venlafaxine_XR$CV_CL),
    Vd_sampled = sample_lognormal_v3(n, p$venlafaxine_XR$Vd, p$venlafaxine_XR$CV_Vd),
    CL_met_sampled = sample_lognormal_v3(n, p$venlafaxine_XR$ODV$CL_odv, p$venlafaxine_XR$ODV$CV_CL_odv),
    Vd_met_sampled = p$venlafaxine_XR$ODV$Vd_odv,
    z_CL = NA_real_, z_Vd = NA_real_, z_CL_met = NA_real_
  )

  bind_rows(sertraline, fluoxetine, paroxetine, ven_ir, ven_xr)
}

# ── 4. 실행 ───────────────────────────────────────────────────────────────────

if (file.exists(RECOVERY_FILE)) {

  cat(sprintf("RECOVERED SAMPLE: %s에서 z값을 역산해 재적용합니다.\n", RECOVERY_FILE))
  sampling_mode <- "recovered_z"

  old_all <- read_csv(RECOVERY_FILE, show_col_types = FALSE)
  verify_v1_consistency(old_all)

  n_per_drug <- old_all %>% count(drug) %>% pull(n) %>% unique()
  if (length(n_per_drug) != 1 || n_per_drug[1] != N) {
    cat(sprintf("참고: 복원 파일의 환자 수(%s)가 N(%d)과 다릅니다 -- 복원 파일 기준 환자 수를 그대로 씁니다.\n",
                paste(n_per_drug, collapse = "/"), N))
  }

  params_all <- bind_rows(
    recover_sertraline_v3(old_all),
    recover_fluoxetine_v3(old_all),
    recover_paroxetine_v3(old_all),
    recover_venlafaxine_v3(old_all, "IR"),
    recover_venlafaxine_v3(old_all, "XR")
  ) %>%
    mutate(sampling_mode = sampling_mode)

} else {

  sampling_mode <- "fresh_seed2025"
  params_all <- fresh_sample_all_v3(N) %>%
    mutate(sampling_mode = sampling_mode)
}

cat(sprintf("Sampling mode: %s\n", sampling_mode))
cat(sprintf("Total rows: %d\n", nrow(params_all)))

# ── 5. D11 재확인용 요약 (본 실행 뒤 7-2 비교에 사용) ────────────────────────
# sertraline 전체, paroxetine NM(=normal)만: 환자별 구현 반감기(0.693*Vd/CL)
# 중앙값이 전형값(27.2h, 21.0h)과 2% 이내여야 함 (HANDOFF_v3 7-2 마지막 문단)

cat("\n-- D11 median-preservation quick check (implemented half-life median) --\n")
params_all %>%
  mutate(t_half_impl = 0.693 * Vd_sampled / CL_sampled) %>%
  filter((drug == "Sertraline") |
         (drug == "Paroxetine" & cyp2d6 == "normal")) %>%
  group_by(drug) %>%
  summarise(t_half_median = round(median(t_half_impl), 2), n = n(), .groups = "drop") %>%
  print()
cat("기대값: Sertraline 27.2h, Paroxetine(normal만) 21.0h -- 2% 이내여야 함 (벗어나면 7-2에서 멈추고 확인)\n")

# ── 6. 저장 ───────────────────────────────────────────────────────────────────

walk(unique(params_all$drug), function(d) {
  df   <- params_all %>% filter(drug == d)
  path <- file.path(OUT_DIR,
                    sprintf("mc_params_v3_%s_n%d.csv",
                            tolower(gsub("_", "", d)), nrow(df)))
  write_csv(df, path)
  cat(sprintf("Saved: %s\n", path))
})

all_path <- file.path(OUT_DIR, sprintf("mc_params_v3_all_n%d.csv", N))
write_csv(params_all, all_path)
cat(sprintf("Saved (all): %s\n", all_path))

cat(sprintf("\nmc_params_v3.R 완료 -- %d rows, sampling_mode=%s\n", nrow(params_all), sampling_mode))
cat("다음 단계: mc_run_v3.R (이 파일의 출력 CSV를 읽어 ODE를 실행)\n\n")
