# =============================================================================
# tables_v3.R
# Phase 6 -- v3 재실행 결과에서 표를 스크립트로 생성 (HANDOFF_TO_NARI_v3.md 7-3)
#
# 만드는 것
#   Table 3 v3 : 지표별(t75/50/25/10%) reached/N 분모를 전부 따로 적음(H4,
#                "t50%의 분모만 적지 마세요")
#   Table 4 v3 : 라벨 parent t1/2 대 모델 t1/2(둘 다 표에 둠, B5) + parent-only
#                약물의 환자별 t50%/구현반감기 중앙값 비 문장
#   Table 5 v3 : 정규화 AUC deficit % (H6: 음수 사례 수를 세서 같이 보고,
#                정의를 바꾸거나 0으로 자르지 않음)
#   Table 6 v3 : paroxetine 표현형별(PM vs non-PM) censoring만 (D5 -- 다른
#                약물은 CYP2D6 모델링 안 함). 거의 없으면 보충표로 옮기라고
#                콘솔에 권고(7-3)
#   Hiemke 2018 비교표 : 환자별 C0가 치료적 참고범위 아래/안/위 드는 비율
#   Michelson 2000 대조표 : 환자별 96h/120h 잔존비율 평균+중앙값 vs Michelson
#   지표별 순위표(H11) : 라벨 t1/2, 모델 t1/2, t50%, 24h reldiff, 절대/정규화
#                AUC deficit 약물 순위 (임상 위험 순위 아님, 각주)
#   Dunn test(H5) : AUC 지표는 전원(모든 환자)으로 재계산 -- v2는 t50%
#                reached 환자만 써서 어긋났음. t50%는 reached만(시간 지표라
#                censored/calc_failed는 값 자체가 없음). Bonferroni로 통일.
#   paroxetine CYP2D6(3.4절) : PM vs non-PM 두 군 Mann-Whitney + rank-biserial r
#                (sensitivity_v3.R과 동일 지표, Kruskal-Wallis 아님 -- 이제
#                두 군뿐이라 RM 없음)
#
# 하지 않는 것(별도 작업 필요)
#   S1(H10, 01_submitted 대비) -- 원 제출본 docx에서 실제 보고된 숫자를
#   읽어와야 함(CSV만으로는 불가). 이 스크립트는 자리만 마련해 둠 --
#   01_submitted/manuscript_SUBMITTED.docx 등에서 숫자를 뽑은 뒤 별도로
#   채울 것.
#
# Input : outputs/mc_v3/mc_results_v3_all_n10000.csv
#         data/pk_parameters/pk_parameters_v3.R
# Output: outputs/tables_v3/ 아래 CSV 여러 개
#
# 실행: source("R/analysis/tables_v3.R")  -- mc_run_v3.R(N=10000) 완료 후
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
})

BASE      <- "~/nari-research/pkpd-antidepressant-sim"
MC3_DIR   <- file.path(BASE, "outputs", "mc_v3")
TABLE_DIR <- file.path(BASE, "outputs", "tables_v3")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "data", "pk_parameters", "pk_parameters_v3.R"))

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

v3_path <- file.path(MC3_DIR, sprintf("mc_results_v3_all_n%d.csv", N))
if (!file.exists(v3_path)) {
  stop(sprintf("%s 없음. mc_params_v3.R -> mc_run_v3.R을 N=%d로 먼저 실행하세요.", v3_path, N))
}
v3 <- read_csv(v3_path, show_col_types = FALSE) %>%
  mutate(drug = factor(drug, levels = drug_levels))

cat(sprintf("v3 rows: %d\n", nrow(v3)))

fmt_mqr <- function(x, digits = 1) {
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf(paste0("%.", digits, "f [%.", digits, "f-%.", digits, "f]"), q[2], q[1], q[3])
}
p_fmt <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.4f", p))

rank_biserial_r <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  U <- suppressWarnings(wilcox.test(x1, x2, exact = FALSE))$statistic
  1 - (2 * U) / (length(x1) * length(x2))
}

# =============================================================================
# TABLE 3 v3 -- 지표별 reached/N 분모를 전부 따로 적음 (H4)
# =============================================================================

metric_denoms <- function(df, status_col, label) {
  tibble(metric = label,
         n_total = nrow(df),
         n_reached = sum(df[[status_col]] == "reached", na.rm = TRUE),
         n_censored = sum(df[[status_col]] == "censored", na.rm = TRUE),
         n_calc_failed = sum(df[[status_col]] == "calc_failed", na.rm = TRUE))
}

table3_denoms <- v3 %>%
  group_by(drug) %>%
  group_modify(~ bind_rows(
    metric_denoms(.x, "status_75pct", "t75%"),
    metric_denoms(.x, "status_50pct", "t50%"),
    metric_denoms(.x, "status_25pct", "t25%"),
    metric_denoms(.x, "status_10pct", "t10%")
  )) %>%
  ungroup() %>%
  mutate(Drug = drug_display[as.character(drug)],
         `reached / N` = sprintf("%d / %d", n_reached, n_total)) %>%
  select(Drug, metric, `reached / N`, n_censored, n_calc_failed)

cat("\n=== Table 3 (분모): reached/N per metric, per drug (H4) ===\n")
print(table3_denoms, n = Inf, width = Inf)
write_csv(table3_denoms, file.path(TABLE_DIR, "table3_v3_denominators.csv"))

table3_metrics <- v3 %>%
  group_by(drug) %>%
  summarise(
    n_total        = n(),
    C0_ngmL        = fmt_mqr(C0 * 1000, 2),
    t75_h          = fmt_mqr(t_75pct[status_75pct == "reached"], 1),
    t50_h          = fmt_mqr(t_50pct[status_50pct == "reached"], 1),
    t25_h          = fmt_mqr(t_25pct[status_25pct == "reached"], 1),
    t10_h          = fmt_mqr(t_10pct[status_10pct == "reached"], 1),
    reldiff_24h    = fmt_mqr(reldiff_24h, 5),
    frac96         = fmt_mqr(frac96, 4),
    frac120        = fmt_mqr(frac120, 4),
    AUCdef_abs     = fmt_mqr(AUC_deficit_7d, 3),
    AUCdef_norm    = fmt_mqr(AUC_deficit_7d_norm_pct, 1),
    median_t50_sort = median(t_50pct[status_50pct == "reached"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median_t50_sort) %>%
  mutate(Rank = row_number(), Drug = drug_display[as.character(drug)]) %>%
  select(Drug, Rank, n_total, `C0 (ng/mL)` = C0_ngmL,
         `t75% (h)` = t75_h, `t50% (h)` = t50_h, `t25% (h)` = t25_h, `t10% (h)` = t10_h,
         `24h reldiff` = reldiff_24h, frac96, frac120,
         `AUC deficit 7d (abs)` = AUCdef_abs, `AUC deficit 7d (norm %)` = AUCdef_norm)

cat("\n=== Table 3: MC PK discontinuation stress metrics (median [Q1-Q3]) ===\n")
print(table3_metrics, n = Inf, width = Inf)
write_csv(table3_metrics, file.path(TABLE_DIR, "table3_v3_mc_metrics.csv"))
cat("각주: 시간 지표(t75/50/25/10%) 중앙값은 reached-only. 분모는 위 table3_v3_denominators.csv 참고.\n")

# =============================================================================
# TABLE 4 v3 -- 라벨 t1/2 대 모델 t1/2 (B5)
# =============================================================================

table4 <- v3 %>%
  group_by(drug) %>%
  summarise(
    t50_median_reached = median(t_50pct[status_50pct == "reached"], na.rm = TRUE),
    t_half_impl_median  = median(t_half_impl_parent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Drug = drug_display[as.character(drug)],
    `Label parent t1/2 (h)` = label_thalf[as.character(drug)],
    `Model implemented t1/2 median (h)` = round(t_half_impl_median, 2),
    `Model t50%C0 median (h)` = round(t50_median_reached, 2),
    `t50%C0 / implemented t1/2` = round(t50_median_reached / t_half_impl_median, 3),
    `Model t1/2 / label t1/2` = round(t_half_impl_median / label_thalf[as.character(drug)], 3)
  ) %>%
  select(Drug, `Label parent t1/2 (h)`, `Model implemented t1/2 median (h)`,
         `Model t1/2 / label t1/2`, `Model t50%C0 median (h)`, `t50%C0 / implemented t1/2`)

cat("\n=== Table 4: label t1/2 vs model t1/2 (both columns, B5) ===\n")
print(table4, n = Inf, width = Inf)
write_csv(table4, file.path(TABLE_DIR, "table4_v3_label_vs_model_thalf.csv"))

parent_only <- table4 %>% filter(Drug %in% c("Sertraline", "Paroxetine"))
cat("\n3.6절에 넣을 문장(parent-only 약물, t50%C0 / 구현반감기 비):\n")
for (i in seq_len(nrow(parent_only))) {
  cat(sprintf("  %s: t50%%C0 median / implemented t1/2 median = %.3f\n",
              parent_only$Drug[i], parent_only$`t50%C0 / implemented t1/2`[i]))
}

# =============================================================================
# TABLE 5 v3 -- 정규화 AUC deficit % (H6: 음수 사례 수를 세서 보고)
# =============================================================================

n_negative <- v3 %>% group_by(drug) %>%
  summarise(n_negative_AUCdef = sum(AUC_deficit_7d_norm_pct < 0, na.rm = TRUE),
            n_total = n(), .groups = "drop")

table5 <- v3 %>%
  group_by(drug) %>%
  summarise(
    C0_ngmL = fmt_mqr(C0 * 1000, 2),
    AUCdef_norm_pct = fmt_mqr(AUC_deficit_7d_norm_pct, 2),
    sort_val = median(AUC_deficit_7d_norm_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(n_negative, by = "drug") %>%
  arrange(desc(sort_val)) %>%
  mutate(Drug = drug_display[as.character(drug)],
         `n negative (H6)` = sprintf("%d / %d", n_negative_AUCdef, n_total)) %>%
  select(Drug, `C0 (ng/mL)` = C0_ngmL, `Normalised AUC deficit 7d (%)` = AUCdef_norm_pct,
         `n negative (H6)`)

cat("\n=== Table 5: normalised AUC deficit (H6: negative count reported, not clipped) ===\n")
print(table5, n = Inf, width = Inf)
write_csv(table5, file.path(TABLE_DIR, "table5_v3_auc_deficit_norm.csv"))

# =============================================================================
# TABLE 6 v3 -- paroxetine 표현형별(PM vs non-PM) censoring만 (D5)
# =============================================================================

table6 <- v3 %>%
  filter(drug == "Paroxetine") %>%
  group_by(cyp2d6) %>%
  summarise(pct_censored_50 = round(mean(status_50pct == "censored") * 100, 3),
            n = n(), .groups = "drop")

cat("\n=== Table 6: paroxetine PM vs non-PM censored % (t50%) ===\n")
print(table6, n = Inf)
write_csv(table6, file.path(TABLE_DIR, "table6_v3_paroxetine_censored_pct.csv"))

max_censor <- max(table6$pct_censored_50, na.rm = TRUE)
if (max_censor < 1) {
  cat(sprintf("\n관찰창 30일에서 censoring이 모두 1%%대 미만(최대 %.3f%%) -- HANDOFF_v3 7-3 권고대로\n", max_censor))
  cat("Table 6을 없애고 본문 한 문장 + 보충표로 옮기는 것을 검토하세요.\n")
} else {
  cat(sprintf("\ncensoring 최대 %.3f%% -- 1%% 넘음, Table 6 유지 권장.\n", max_censor))
}

# =============================================================================
# Hiemke 2018 (AGNP) 치료적 참고범위 비교 -- 환자별 C0(최저농도) 분포
# =============================================================================

hiemke_ranges <- tribble(
  ~drug,           ~lo,  ~hi,  ~note,
  "Sertraline",     10,  150,  "parent만",
  "Paroxetine",     20,   65,  "parent만",
  "Fluoxetine",    120,  500,  "parent+norfluoxetine 합",
  "Venlafaxine_IR",100,  400,  "parent+ODV 합",
  "Venlafaxine_XR",100,  400,  "parent+ODV 합"
)

# 합산 지표가 필요한 fluoxetine/venlafaxine은 mc_run_v3.R의 conc_col이 이미
# Ctotal(=Cp+Cm) 기준이므로 C0 컬럼 자체가 이미 "합"임 -- 별도 합산 불필요.
hiemke_table <- v3 %>%
  left_join(hiemke_ranges, by = "drug") %>%
  filter(!is.na(lo)) %>%
  group_by(drug) %>%
  summarise(
    n = n(),
    pct_below = round(mean(C0 * 1000 < lo) * 100, 2),
    pct_within = round(mean(C0 * 1000 >= lo & C0 * 1000 <= hi) * 100, 2),
    pct_above = round(mean(C0 * 1000 > hi) * 100, 2),
    range_lo = first(lo), range_hi = first(hi), note = first(note),
    .groups = "drop"
  ) %>%
  mutate(Drug = drug_display[as.character(drug)],
         `Hiemke range (ng/mL)` = sprintf("%d-%d (%s)", range_lo, range_hi, note)) %>%
  select(Drug, `Hiemke range (ng/mL)`, n,
         `% below` = pct_below, `% within` = pct_within, `% above` = pct_above)

# paroxetine PM은 범위 위로 가는 것이 예상(HANDOFF_v3 7-3) -- 별도로도 보고
hiemke_paroxetine_by_phenotype <- v3 %>%
  filter(drug == "Paroxetine") %>%
  left_join(hiemke_ranges, by = "drug") %>%
  group_by(cyp2d6) %>%
  summarise(n = n(),
            pct_below = round(mean(C0 * 1000 < lo) * 100, 2),
            pct_within = round(mean(C0 * 1000 >= lo & C0 * 1000 <= hi) * 100, 2),
            pct_above = round(mean(C0 * 1000 > hi) * 100, 2),
            .groups = "drop")

cat("\n=== Hiemke 2018 (AGNP) therapeutic reference range comparison ===\n")
cat("(치료적 참고범위이며 특정 용량의 예상 농도 구간이 아님 -- 판정 아니고 비율만 보고)\n")
print(hiemke_table, n = Inf, width = Inf)
cat("\n-- paroxetine, phenotype별 (PM은 범위 위로 가는 것이 예상) --\n")
print(hiemke_paroxetine_by_phenotype, n = Inf)
write_csv(hiemke_table, file.path(TABLE_DIR, "hiemke_v3_reference_range.csv"))
write_csv(hiemke_paroxetine_by_phenotype, file.path(TABLE_DIR, "hiemke_v3_paroxetine_by_phenotype.csv"))

# =============================================================================
# Michelson 2000 대조 -- 환자별 96h/120h 잔존비율 평균+중앙값
# =============================================================================

michelson_ref <- tribble(
  ~drug_group,                              ~michelson_remaining_frac,
  "Fluoxetine + norfluoxetine",              0.703,
  "Sertraline + desmethylsertraline",        0.265,
  "Paroxetine",                              0.133
)

michelson_table <- v3 %>%
  group_by(drug) %>%
  summarise(
    n = n(),
    mean_frac96   = round(mean(frac96, na.rm = TRUE), 4),
    median_frac96 = round(median(frac96, na.rm = TRUE), 4),
    mean_frac120   = round(mean(frac120, na.rm = TRUE), 4),
    median_frac120 = round(median(frac120, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  mutate(Drug = drug_display[as.character(drug)])

cat("\n=== Michelson 2000 comparison (96h/120h remaining fraction, model vs observed) ===\n")
cat("(Michelson: 5일 위약대체 후 환자별 감소율 평균. sertraline은 모델이 parent만이라\n")
cat(" 직접 비교 불가 -- 대사체 합 관찰값과 다름. paroxetine은 PM 혼입 전체값.)\n")
print(michelson_table, n = Inf, width = Inf)
print(michelson_ref)
write_csv(michelson_table, file.path(TABLE_DIR, "michelson_v3_model.csv"))
write_csv(michelson_ref, file.path(TABLE_DIR, "michelson_v3_reference.csv"))

# =============================================================================
# 지표별 순위표 (H11) -- 임상 위험 순위 아님
# =============================================================================

rank_table <- v3 %>%
  group_by(drug) %>%
  summarise(
    model_thalf    = median(t_half_impl_parent, na.rm = TRUE),
    t50_median     = median(t_50pct[status_50pct == "reached"], na.rm = TRUE),
    reldiff24_median = median(reldiff_24h, na.rm = TRUE),
    AUCdef_abs_median  = median(AUC_deficit_7d, na.rm = TRUE),
    AUCdef_norm_median = median(AUC_deficit_7d_norm_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Drug = drug_display[as.character(drug)],
    label_thalf_val = label_thalf[as.character(drug)],
    rank_label_thalf   = rank(label_thalf_val),
    rank_model_thalf   = rank(model_thalf),
    rank_t50           = rank(t50_median),
    rank_reldiff24      = rank(-abs(reldiff24_median)),   # 빠르게 떨어질수록 순위 1
    rank_AUCdef_abs     = rank(-AUCdef_abs_median),
    rank_AUCdef_norm    = rank(-AUCdef_norm_median)
  ) %>%
  select(Drug, label_thalf_val, rank_label_thalf, model_thalf, rank_model_thalf,
         t50_median, rank_t50, reldiff24_median, rank_reldiff24,
         AUCdef_abs_median, rank_AUCdef_abs, AUCdef_norm_median, rank_AUCdef_norm)

cat("\n=== Rank table (H11) -- NOT a clinical risk ranking (footnote) ===\n")
print(rank_table, n = Inf, width = Inf)
write_csv(rank_table, file.path(TABLE_DIR, "rank_table_v3.csv"))

# =============================================================================
# Dunn test (H5) -- AUC 지표는 전원(모든 환자)으로 재계산. t50%는 reached만.
# =============================================================================

dunn_t50_input <- v3 %>% filter(status_50pct == "reached") %>%
  mutate(drug = factor(drug, levels = drug_levels))
dunn_auc_input <- v3 %>% filter(!is.na(AUC_deficit_7d)) %>%        # H5: 전원(reached 필터 없음)
  mutate(drug = factor(drug, levels = drug_levels))

cat(sprintf("\nDunn test (t50%%): n=%d (reached only, 시간 지표라 불가피)\n", nrow(dunn_t50_input)))
cat(sprintf("Dunn test (AUC):  n=%d (전원, H5 수정 -- v2는 reached만 써서 Methods와 어긋났음)\n", nrow(dunn_auc_input)))

dunn_t50 <- dunn_t50_input %>%
  dunn_test(t_50pct ~ drug, p.adjust.method = "bonferroni") %>%
  select(group1, group2, statistic, p, p.adj, p.adj.signif) %>%
  mutate(metric = "t50% C0 (h, reached only)", statistic = round(statistic, 3),
         p = p_fmt(p), p.adj = p_fmt(p.adj))

dunn_auc <- dunn_auc_input %>%
  dunn_test(AUC_deficit_7d ~ drug, p.adjust.method = "bonferroni") %>%
  select(group1, group2, statistic, p, p.adj, p.adj.signif) %>%
  mutate(metric = "AUC deficit 7d (abs, all patients)", statistic = round(statistic, 3),
         p = p_fmt(p), p.adj = p_fmt(p.adj))

dunn_auc_norm <- dunn_auc_input %>%
  dunn_test(AUC_deficit_7d_norm_pct ~ drug, p.adjust.method = "bonferroni") %>%
  select(group1, group2, statistic, p, p.adj, p.adj.signif) %>%
  mutate(metric = "AUC deficit 7d (norm %, all patients)", statistic = round(statistic, 3),
         p = p_fmt(p), p.adj = p_fmt(p.adj))

dunn_all <- bind_rows(dunn_t50, dunn_auc, dunn_auc_norm) %>% arrange(metric, group1, group2)

cat("\n=== Dunn test (S2, Bonferroni-adjusted) ===\n")
print(dunn_all, n = Inf, width = Inf)
write_csv(dunn_all, file.path(TABLE_DIR, "s2_dunn_v3_bonferroni.csv"))

# =============================================================================
# paroxetine CYP2D6 (3.4절) -- PM vs non-PM, Mann-Whitney + rank-biserial r
# =============================================================================

par_t50 <- v3 %>% filter(drug == "Paroxetine", status_50pct == "reached")

n_pm  <- sum(par_t50$cyp2d6 == "poor")
n_nm  <- sum(par_t50$cyp2d6 == "normal")

wtest <- suppressWarnings(wilcox.test(
  t_50pct ~ cyp2d6,
  data = par_t50 %>% filter(cyp2d6 %in% c("poor", "normal")),
  exact = FALSE
))
eff_r <- rank_biserial_r(par_t50$t_50pct[par_t50$cyp2d6 == "poor"],
                          par_t50$t_50pct[par_t50$cyp2d6 == "normal"])

cyp_summary <- par_t50 %>%
  group_by(cyp2d6) %>%
  summarise(n = n(), median_t50_h = round(median(t_50pct), 2), .groups = "drop")

cat("\n=== 3.4절: paroxetine CYP2D6 (PM vs non-PM, two groups only -- D5) ===\n")
cat(sprintf("n: PM=%d, non-PM(normal)=%d (실제 n, 'n=10,000 per group' 식 표현 금지, H5/3.4절)\n", n_pm, n_nm))
print(cyp_summary)
cat(sprintf("median ratio PM/non-PM = %.3f\n",
            cyp_summary$median_t50_h[cyp_summary$cyp2d6 == "poor"] /
            cyp_summary$median_t50_h[cyp_summary$cyp2d6 == "normal"]))
cat(sprintf("Mann-Whitney U = %.1f, p = %s, rank-biserial r = %.3f\n",
            wtest$statistic, p_fmt(wtest$p.value), eff_r))

cyp_out <- cyp_summary %>%
  mutate(mann_whitney_p = p_fmt(wtest$p.value), rank_biserial_r = round(eff_r, 3))
write_csv(cyp_out, file.path(TABLE_DIR, "section3_4_v3_paroxetine_cyp2d6.csv"))

# =============================================================================
# S1 (H10) -- 추가 지시 1(2-4): 원 제출본 실제 수치(O_SUBMISSION_NUMBERS.md,
#   덕현 쪽 Claude가 01_submitted/manuscript_no_titlepage.docx와
#   supplementary.docx 표 XML에서 직접 추출) 대 v3 결과.
#   주의: 원 제출본 첫 열은 "Css"(투여간격 평균, 중단 전)이고 v3의 C0(중단
#   시점 농도)와 정의가 다름 -- 직접 비교하지 말고 각자 표기.
#   원 제출본 Table 2는 t50% >= 1h 제외 후 값이며, 제외 후 n이 약물마다
#   다름(보충자료 표). v3에도 같은 기준(status_50pct=="reached" & t_50pct>=1)
#   을 적용해 "같은 정의"로 만든 뒤 비교.
# =============================================================================

cat("\n=== S1 (H10): 원 제출본 대 v3 비교 ===\n")

# 원 제출본 숫자 (O_SUBMISSION_NUMBERS.md, 01_submitted docx에서 직접 추출,
# 덕현 쪽 Claude 작업, 2026-10-28)
s1_submitted <- tribble(
  ~drug,              ~rank_submitted, ~Css_ngmL, ~t50_h_submitted, ~t25_h_submitted, ~reldiff24_submitted, ~AUCdef_submitted, ~n_main, ~n_noexclusion,
  "Venlafaxine_XR",   1,  72.5,  7.2,  21.5, 0.0324, 11.09, 7417, 8004,
  "Venlafaxine_IR",   2,  65.0, 10.3,  22.8, 0.0319,  9.84, 9901, 9970,
  "Paroxetine",       3,  27.4, 16.4,  43.0, 0.0245,  3.71, 7918, 8195,
  "Sertraline",       4,  61.1, 17.1,  44.4, 0.0242,  8.17, 9373, 9524,
  "Fluoxetine",       5, 137.8, 349,   649,  0.0026,  3.59, 9710, 9710
)

# v3 쪽: 원 제출본과 같은 "t50% >= 1h, 도달(reached)" 기준 적용
s1_v3 <- v3 %>%
  filter(status_50pct == "reached", t_50pct >= 1) %>%
  group_by(drug) %>%
  summarise(
    n_v3             = n(),
    C0_v3            = median(C0, na.rm = TRUE) * 1000,   # ng/mL
    t50_h_v3         = median(t_50pct, na.rm = TRUE),
    t25_h_v3         = median(t_25pct[status_25pct == "reached"], na.rm = TRUE),
    reldiff24_v3     = median(reldiff_24h, na.rm = TRUE),
    AUCdef_v3        = median(AUC_deficit_7d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(rank_v3 = rank(t50_h_v3))

s1_table <- s1_submitted %>%
  left_join(s1_v3, by = "drug") %>%
  mutate(
    n_excluded_pct = round(100 * (1 - n_main / n_noexclusion), 2),
    note = "Css(submitted, pre-discontinuation interval-average) vs C0(v3, at-discontinuation) are NOT the same definition -- do not difference them directly"
  ) %>%
  arrange(rank_submitted)

cat("원 제출본 대 v3 (둘 다 t50%>=1h 제외 기준):\n")
print(s1_table, n = Inf, width = Inf)

s1_rho <- suppressWarnings(cor.test(s1_table$rank_submitted, s1_table$rank_v3, method = "spearman"))
cat(sprintf("\nSpearman rho (원 제출본 순위 vs v3 순위, t50%% 기준) = %.3f (S=%.1f, p=%s)\n",
            s1_rho$estimate, s1_rho$statistic, p_fmt(s1_rho$p.value)))
cat("정성적 과장 표현('순위가 전혀 바뀌지 않았다' 등) 대신 이 rho 값을 그대로 보고할 것\n")
cat("(HANDOFF_v3 7-3 S1 항목, 추가 지시 1 2-4).\n\n")

write_csv(s1_table, file.path(TABLE_DIR, "S1_v3_vs_submitted.csv"))

cat("tables_v3.R 완료. outputs/tables_v3/ 아래 CSV 다수 생성.\n")
cat("다음 단계: Hiemke/Michelson/rank/Dunn 결과를 들고 7-4(원고 수정)로.\n\n")
