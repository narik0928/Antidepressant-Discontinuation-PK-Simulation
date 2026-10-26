# =============================================================================
# make_tables_v2.R
# Phase 6 (RESUBMISSION_PLAN.md 3단계-4) — v2 재실행 결과에서 표를 스크립트로
# 생성. 숫자를 손으로 옮기지 않기 위함.
#
# 만드는 것 (RESUBMISSION_PLAN.md 3단계-4 요구사항):
#   - Table 3 v2 : MC PK discontinuation stress metrics (median [Q1–Q3], 순위)
#   - Table 4 v2 : 시뮬레이션 t50%C0 vs 라벨 parent t½ 비율
#   - Table 5 v2 : 정규화 AUC deficit % (A2 버그 수정 반영)
#   - Table 6 v2 : CYP2D6 phenotype별 censored 비율 (v1의 "<1h 제외 비율" 대체 —
#                  <1h 제외 규칙 자체가 v2에서 없어졌으므로 의미가 바뀜)
#   - S1         : 구 정의(Css, v1) vs 신 정의(C0, v2) 환자별 짝지은 비교
#   - B1         : 5개 regimen의 t50%C0 순위 비교 (v1 vs v2)
#   - Dunn test  : t50%C0 / AUC_deficit_7d 약물 간 비교, raw + Bonferroni adjusted p
#                  (phase5_statistics.R과 동일한 방법론, 단 <1h 제외는 없앰 —
#                  v2는 애초에 그런 환자가 거의 없음, censored만 별도 처리)
#
# Input : outputs/mc_v2/mc_results_v2_all_n10000.csv
#         outputs/mc/mc_results_all_n10000.csv (v1, S1/B1 비교용)
#         outputs/qualification_v2/model_qualification_table_v2.csv
#         outputs/qualification_v2/parent_metabolite_ratio_v2.csv
#         data/pk_parameters/pk_parameters_v2.R (라벨 t½ 값)
# Output: outputs/tables_v2/ 아래 CSV 7개 + 콘솔에 사람이 읽기 좋은 형태로 출력
#
# 실행: source("R/analysis/make_tables_v2.R")
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
})

BASE      <- "~/nari-research/pkpd-antidepressant-sim"
MC_DIR    <- file.path(BASE, "outputs", "mc")
MC2_DIR   <- file.path(BASE, "outputs", "mc_v2")
QUAL2_DIR <- file.path(BASE, "outputs", "qualification_v2")
TABLE_DIR <- file.path(BASE, "outputs", "tables_v2")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "data", "pk_parameters", "pk_parameters_v2.R"))

N <- 10000

drug_levels  <- c("Sertraline", "Fluoxetine", "Paroxetine", "Venlafaxine_IR", "Venlafaxine_XR")
drug_display <- c(Sertraline = "Sertraline", Fluoxetine = "Fluoxetine + NFX",
                   Paroxetine = "Paroxetine", Venlafaxine_IR = "Venlafaxine IR",
                   Venlafaxine_XR = "Venlafaxine XR")
label_thalf  <- c(Sertraline = pk_params$sertraline$t_half,
                   Fluoxetine = pk_params$fluoxetine$t_half,
                   Paroxetine = pk_params$paroxetine$t_half,
                   Venlafaxine_IR = pk_params$venlafaxine_IR$t_half,
                   Venlafaxine_XR = pk_params$venlafaxine_XR$t_half)

# ── 0. Load ───────────────────────────────────────────────────────────────────
v2 <- read_csv(file.path(MC2_DIR, sprintf("mc_results_v2_all_n%d.csv", N)), show_col_types = FALSE) %>%
  mutate(drug = factor(drug, levels = drug_levels))
v1 <- read_csv(file.path(MC_DIR, sprintf("mc_results_all_n%d.csv", N)), show_col_types = FALSE) %>%
  mutate(drug = factor(drug, levels = drug_levels))

cat(sprintf("v2 rows: %d | v1 rows: %d\n", nrow(v2), nrow(v1)))

# ── helpers ───────────────────────────────────────────────────────────────────
fmt_mqr <- function(x, digits = 1) {
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf(paste0("%.", digits, "f [%.", digits, "f–%.", digits, "f]"), q[2], q[1], q[3])
}

p_fmt <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.4f", p))

# ============================================================================
# TABLE 3 v2 — MC PK discontinuation stress metrics
# v1과 달리 <1h 제외 없음. 시간 기반 지표(t50%, t25%)는 reached 환자만 요약에
# 넣고, n_reached/N을 같이 적어서 실제 분모를 명시 (RESUBMISSION_PLAN "Table
# 6(실제 분모)" 요구사항을 Table 3에도 동일하게 적용 — 투명성을 위해).
# pct_drop/reldiff/AUC deficit는 모든 환자(C0>0)에서 계산 가능하므로 전원 포함.
# ============================================================================

table3 <- v2 %>%
  group_by(drug) %>%
  summarise(
    n_total       = n(),
    n_reached_t50 = sum(status_50pct == "reached"),
    C0_ngmL       = fmt_mqr(C0 * 1000, 1),
    t50_h         = fmt_mqr(t_50pct[status_50pct == "reached"], 1),
    t25_h         = fmt_mqr(t_25pct[status_25pct == "reached"], 1),
    abs_reldiff_24h = fmt_mqr(abs(reldiff_24h), 4),
    AUCdef_7d_mgLh  = fmt_mqr(AUC_deficit_7d, 2),
    median_t50_sort = median(t_50pct[status_50pct == "reached"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median_t50_sort) %>%
  mutate(Rank = row_number(),
         Drug = drug_display[as.character(drug)],
         `n reached / N` = sprintf("%d / %d", n_reached_t50, n_total)) %>%
  select(Drug, Rank, `n reached / N`, `C0 (ng/mL)` = C0_ngmL,
         `t50% C0 (h)` = t50_h, `t25% C0 (h)` = t25_h,
         `|Rel. dC/dt| at 24h (h-1)` = abs_reldiff_24h,
         `AUC deficit 7d (mg/L*h)` = AUCdef_7d_mgLh)

cat("\n========== Table 3 v2: MC PK discontinuation stress metrics ==========\n")
print(table3, n = Inf, width = Inf)
write_csv(table3, file.path(TABLE_DIR, "table3_v2_mc_metrics.csv"))

# ============================================================================
# TABLE 4 v2 — Simulated t50%C0 vs label parent t1/2
# ============================================================================

table4 <- v2 %>%
  filter(status_50pct == "reached") %>%
  group_by(drug) %>%
  summarise(t50_median = median(t_50pct, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Drug = drug_display[as.character(drug)],
    `Label parent t1/2 (h)` = label_thalf[as.character(drug)],
    `Simulated median t50% C0 (h)` = round(t50_median, 1),
    `Ratio (simulated / label)` = round(t50_median / label_thalf[as.character(drug)], 2)
  ) %>%
  select(Drug, `Label parent t1/2 (h)`, `Simulated median t50% C0 (h)`, `Ratio (simulated / label)`)

cat("\n========== Table 4 v2: simulated t50%C0 vs label t1/2 ==========\n")
print(table4, n = Inf, width = Inf)
write_csv(table4, file.path(TABLE_DIR, "table4_v2_label_ratio.csv"))

# ============================================================================
# TABLE 5 v2 — Normalised AUC deficit (%) — A2 fix: 이제 환자별로 정규화된
# 진짜 % 값을 바로 median [Q1-Q3]로 요약 (v1처럼 drug-level Css/7일로 역산
# 하지 않음 — 환자마다 C0가 다르므로 patient-level 정규화가 맞는 방식).
# ============================================================================

table5 <- v2 %>%
  group_by(drug) %>%
  summarise(
    C0_ngmL   = fmt_mqr(C0 * 1000, 1),
    AUCdef_norm_pct = fmt_mqr(AUC_deficit_7d_norm_pct, 1),
    sort_val = median(AUC_deficit_7d_norm_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sort_val)) %>%
  mutate(Drug = drug_display[as.character(drug)]) %>%
  select(Drug, `C0 (ng/mL)` = C0_ngmL, `Normalised AUC deficit 7d (%)` = AUCdef_norm_pct)

cat("\n========== Table 5 v2: normalised AUC deficit (A2 fix applied) ==========\n")
print(table5, n = Inf, width = Inf)
write_csv(table5, file.path(TABLE_DIR, "table5_v2_auc_deficit_norm.csv"))

# ============================================================================
# TABLE 6 v2 — CYP2D6 phenotype-stratified CENSORED percentage
# v1의 "<1h 제외 비율" 표는 그 제외 규칙 자체가 사라져서 더 이상 성립하지
# 않음. v2는 대신 "관찰 종료까지 50% 미도달(censored) 비율"을 같은 형식으로
# 보고 — 숫자가 작게 나올수록 D1/D3/D4 수정이 관찰창 문제를 해결했다는 뜻.
# ============================================================================

table6 <- v2 %>%
  filter(drug != "Sertraline") %>%
  group_by(drug, cyp2d6) %>%
  summarise(pct_censored = mean(status_50pct == "censored") * 100, .groups = "drop") %>%
  pivot_wider(names_from = cyp2d6, values_from = pct_censored, values_fill = 0)

table6_overall <- v2 %>%
  filter(drug != "Sertraline") %>%
  group_by(drug) %>%
  summarise(Overall = mean(status_50pct == "censored") * 100, .groups = "drop")

sert_row <- v2 %>% filter(drug == "Sertraline") %>%
  summarise(Overall = mean(status_50pct == "censored") * 100)

table6_final <- table6 %>%
  left_join(table6_overall, by = "drug") %>%
  mutate(Drug = drug_display[as.character(drug)]) %>%
  select(Drug, poor, normal, rapid, Overall) %>%
  bind_rows(tibble(
    Drug = "Sertraline (no CYP2D6 stratification modelled)",
    poor = NA, normal = NA, rapid = NA, Overall = sert_row$Overall
  )) %>%
  mutate(across(c(poor, normal, rapid, Overall), ~ ifelse(is.na(.x), "—", sprintf("%.2f%%", .x))))

cat("\n========== Table 6 v2: CYP2D6-stratified CENSORED % (not '<1h excluded' -- that rule is gone) ==========\n")
print(table6_final, n = Inf, width = Inf)
write_csv(table6_final, file.path(TABLE_DIR, "table6_v2_censored_pct.csv"))

# ============================================================================
# S1 — 구 정의(Css, v1) vs 신 정의(C0, v2) 환자별 짝지은 비교
# 같은 코호트를 재사용했으므로 (drug, id)로 정확히 짝지을 수 있음.
#
# 중요: v1→v2는 정의(D1: Css→C0)와 F(D3: F<1→F=1)가 동시에 바뀌었으므로,
# C0(v2) vs Css(v1)만 비교하면 두 효과가 섞여서 "정의를 바꿔서 생긴 차이"라고
# 말할 수 없음 (RESUBMISSION_PLAN.md 0절: "같은 가상 환자끼리 짝지어 비교해야
# 원인을 나눌 수 있다"). 그래서 세 가지를 분리해서 보고:
#   (1) 순수 D1 효과  : C0(v2) vs legacy_Css(v2) -- 둘 다 F=1, 정의만 다름
#   (2) 순수 D3 효과  : legacy_Css(v2) vs Css(v1) -- 정의(평균 방식)는 같고 F만 다름
#   (3) 합산 효과     : C0(v2) vs Css(v1) -- 실제 원고에 들어갈 숫자가 바뀐 폭
# ============================================================================

paired <- v1 %>%
  select(drug, id, CL_sampled, Vd_sampled, Css_v1 = Css, t50_v1 = t_50pct, AUCdef_v1 = AUC_deficit_7d) %>%
  inner_join(
    v2 %>% select(drug, id, C0_v2 = C0, legacy_Css_v2 = legacy_Css,
                  t50_v2 = t_50pct, status50_v2 = status_50pct, AUCdef_v2 = AUC_deficit_7d),
    by = c("drug", "id")
  )

cat(sprintf("\nS1 짝지은 환자 수: %d (v1과 v2 모두 존재하는 (drug,id) 조합)\n", nrow(paired)))

s1 <- paired %>%
  group_by(drug) %>%
  summarise(
    n_paired              = n(),
    Css_v1_ngmL_median    = round(median(Css_v1, na.rm = TRUE) * 1000, 1),
    legacy_Css_v2_ngmL_median = round(median(legacy_Css_v2, na.rm = TRUE) * 1000, 1),
    C0_v2_ngmL_median     = round(median(C0_v2, na.rm = TRUE) * 1000, 1),
    # (1) 순수 D1 효과 (같은 F=1 안에서 정의만 바뀜)
    pct_C0_lower_than_legacyCss_sameF = round(mean(C0_v2 < legacy_Css_v2, na.rm = TRUE) * 100, 1),
    # (2) 순수 D3 효과 (같은 평균-방식 정의 안에서 F만 바뀜)
    pct_legacyCss_higher_than_Css_v1 = round(mean(legacy_Css_v2 > Css_v1, na.rm = TRUE) * 100, 1),
    # (3) 합산 효과 (실제 원고 숫자가 움직인 폭)
    pct_C0_lower_than_Css_v1 = round(mean(C0_v2 < Css_v1, na.rm = TRUE) * 100, 1),
    t50_v1_median         = round(median(t50_v1, na.rm = TRUE), 1),
    t50_v2_median_reached = round(median(t50_v2[status50_v2 == "reached"], na.rm = TRUE), 1),
    n_reached_v2          = sum(status50_v2 == "reached"),
    n_censored_v2         = sum(status50_v2 == "censored"),
    n_NA_v1               = sum(is.na(t50_v1)),
    spearman_rho_t50      = suppressWarnings(cor(t50_v1, t50_v2, method = "spearman", use = "complete.obs")),
    .groups = "drop"
  ) %>%
  mutate(Drug = drug_display[as.character(drug)]) %>%
  relocate(Drug)

cat("\n========== S1: old (Css, v1) vs new (C0, v2) definition, paired by patient ==========\n")
cat("(D1 효과와 D3 효과를 분리해서 보고 -- 합산 효과만 보면 두 변경이 섞여 보임)\n")
print(s1, n = Inf, width = Inf)
write_csv(s1, file.path(TABLE_DIR, "s1_v2_old_vs_new_paired.csv"))
write_csv(paired, file.path(TABLE_DIR, "s1_v2_old_vs_new_paired_patientlevel.csv"))

# ============================================================================
# B1 — 5개 regimen 순위 비교: v1 vs v2 (t50%C0 median 기준, 빠른 순)
# ============================================================================

rank_v1 <- v1 %>%
  group_by(drug) %>%
  summarise(t50_median_v1 = median(t_50pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(t50_median_v1) %>%
  mutate(rank_v1 = row_number())

rank_v2 <- v2 %>%
  filter(status_50pct == "reached") %>%
  group_by(drug) %>%
  summarise(t50_median_v2 = median(t_50pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(t50_median_v2) %>%
  mutate(rank_v2 = row_number())

b1 <- rank_v1 %>%
  inner_join(rank_v2, by = "drug") %>%
  mutate(
    Drug = drug_display[as.character(drug)],
    rank_changed = rank_v1 != rank_v2,
    t50_median_v1 = round(t50_median_v1, 1),
    t50_median_v2 = round(t50_median_v2, 1)
  ) %>%
  arrange(rank_v2) %>%
  select(Drug, rank_v1, t50_median_v1, rank_v2, t50_median_v2, rank_changed)

cat("\n========== B1: regimen ranking comparison (v1 vs v2, by t50% C0 median) ==========\n")
print(b1, n = Inf, width = Inf)
write_csv(b1, file.path(TABLE_DIR, "b1_v2_ranking_comparison.csv"))

if (any(b1$rank_changed)) {
  cat("\n*** 순위가 바뀐 regimen이 있음 -- Discussion의 순위 서술을 반드시 다시 확인 ***\n")
} else {
  cat("\n순위 자체는 v1과 v2가 동일함 (절대값은 달라짐).\n")
}

# ============================================================================
# Dunn test (t50%C0, AUC_deficit_7d) — 약물 간 비교, v2 기준
# v1(phase5_statistics.R)과 달리 <1h 제외 없음. t50%는 reached 환자만 포함
# (censored/calc_failed는 정의상 구간 비교가 불가능하므로 제외하되, 그 수가
# 얼마나 되는지는 위 Table 3/6에 이미 보고함 -- 숨기지 않음).
# ============================================================================

dunn_input <- v2 %>%
  filter(status_50pct == "reached", !is.na(AUC_deficit_7d)) %>%
  mutate(drug = factor(drug, levels = drug_levels))

cat(sprintf("\nDunn test 대상 환자 수: %d (reached, AUC_deficit_7d non-NA)\n", nrow(dunn_input)))

dunn_t50 <- dunn_input %>%
  dunn_test(t_50pct ~ drug, p.adjust.method = "bonferroni") %>%
  select(group1, group2, statistic, p, p.adj, p.adj.signif) %>%
  mutate(metric = "t50% C0 (h)", statistic = round(statistic, 3),
         p = p_fmt(p), p.adj = p_fmt(p.adj))

dunn_auc <- dunn_input %>%
  dunn_test(AUC_deficit_7d ~ drug, p.adjust.method = "bonferroni") %>%
  select(group1, group2, statistic, p, p.adj, p.adj.signif) %>%
  mutate(metric = "AUC deficit 7d (abs)", statistic = round(statistic, 3),
         p = p_fmt(p), p.adj = p_fmt(p.adj))

dunn_auc_norm <- dunn_input %>%
  dunn_test(AUC_deficit_7d_norm_pct ~ drug, p.adjust.method = "bonferroni") %>%
  select(group1, group2, statistic, p, p.adj, p.adj.signif) %>%
  mutate(metric = "AUC deficit 7d (normalised %)", statistic = round(statistic, 3),
         p = p_fmt(p), p.adj = p_fmt(p.adj))

dunn_all <- bind_rows(dunn_t50, dunn_auc, dunn_auc_norm) %>%
  arrange(metric, group1, group2)

cat("\n========== Dunn test (Bonferroni-adjusted), v2 ==========\n")
print(dunn_all, n = Inf, width = Inf)
write_csv(dunn_all, file.path(TABLE_DIR, "dunn_v2_drug_comparison.csv"))

# CYP2D6 Kruskal-Wallis effect (동일 방법론, v2 데이터)
kw_cyp <- dunn_input %>%
  filter(drug != "Sertraline") %>%
  group_by(drug) %>%
  do(
    kw      = kruskal_test(., t_50pct ~ cyp2d6),
    effsize = kruskal_effsize(., t_50pct ~ cyp2d6)
  ) %>%
  summarise(
    drug      = first(drug),
    H         = round(kw$statistic, 2),
    df        = kw$df,
    p_value   = p_fmt(kw$p),
    eta_sq    = round(effsize$effsize, 3),
    magnitude = effsize$magnitude,
    .groups   = "drop"
  )

cat("\n========== Kruskal-Wallis: CYP2D6 effect on t50% C0, v2 ==========\n")
print(kw_cyp, n = Inf, width = Inf)
write_csv(kw_cyp, file.path(TABLE_DIR, "kruskal_v2_cyp2d6.csv"))

# ============================================================================
# 요약
# ============================================================================
cat("\n✅ make_tables_v2.R 완료. outputs/tables_v2/ 아래 9개 CSV 생성:\n")
cat("  table3_v2_mc_metrics.csv, table4_v2_label_ratio.csv, table5_v2_auc_deficit_norm.csv,\n")
cat("  table6_v2_censored_pct.csv, s1_v2_old_vs_new_paired.csv (+_patientlevel.csv),\n")
cat("  b1_v2_ranking_comparison.csv, dunn_v2_drug_comparison.csv, kruskal_v2_cyp2d6.csv\n")
cat("\n이 파일들을 나한테 보내주면 원고 Table 3-6/S1/B1을 이 숫자로 다시 쓸게.\n")
