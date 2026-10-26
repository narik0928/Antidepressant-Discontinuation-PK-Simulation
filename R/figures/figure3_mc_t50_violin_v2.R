# =============================================================================
# figure3_mc_t50_violin_v2.R
# Figure 3 (v2): Monte Carlo — t50% C0 distribution by drug and CYP2D6 phenotype
#
# Input : outputs/mc_v2/mc_results_v2_all_n10000.csv
# Output: outputs/figures/figure3_mc_t50_violin_v2.png
#
# v2 changes: D1/D3/D4 — reads v2 MC output; no t_50pct >= 1 exclusion
# (status_50pct distinguishes reached/censored instead, see Table 3/6);
# fluoxetine's extended 400-day observation window pushes reached values
# up to ~9,500 h, so the y-axis range/breaks are widened accordingly.
#
# Design:
#   Violin (distribution) + Boxplot (IQR) per drug × CYP2D6 phenotype
#   Log10 y-axis — Fluoxetine t50% median ~368h vs Venlafaxine XR ~13h
#   Drug order: ascending median t50% (fastest decline left)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ── 0. Setup ──────────────────────────────────────────────────────────────────
BASE    <- "~/nari-research/pkpd-antidepressant-sim"
FIG_DIR <- file.path(BASE, "outputs", "figures")
OUT_DIR <- file.path(BASE, "outputs", "mc_v2")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

N <- 10000

# ── 1. Load results ───────────────────────────────────────────────────────────
results <- read_csv(
  file.path(OUT_DIR, sprintf("mc_results_v2_all_n%d.csv", N)),
  show_col_types = FALSE
)

# ── 2. Status report (reached / censored / calc_failed) ───────────────────────
status_report <- results %>%
  group_by(drug, status_50pct) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(drug) %>%
  mutate(pct = round(n / sum(n) * 100, 2)) %>%
  ungroup()
cat("\n── t50% status report (reached / censored / calc_failed) ──\n")
print(status_report, n = Inf)

# ── 3. Data prep ──────────────────────────────────────────────────────────────
# Drug order: ascending median t50% (fastest elimination first)
drug_order <- c("Venlafaxine_XR", "Venlafaxine_IR",
                "Paroxetine", "Sertraline", "Fluoxetine")

drug_labels <- c(
  "Venlafaxine_XR" = "Venlafaxine\nXR",
  "Venlafaxine_IR" = "Venlafaxine\nIR",
  "Paroxetine"     = "Paroxetine",
  "Sertraline"     = "Sertraline\n(CYP2D6 N/A)",
  "Fluoxetine"     = "Fluoxetine\n+NFX"
)

cyp2d6_colors <- c(
  "poor"   = "#C0392B",   # 빨강 — CL 낮음, 농도 높음, 소실 느림
  "normal" = "#2980B9",   # 파랑
  "rapid"  = "#27AE60"    # 초록 — CL 높음, 농도 낮음, 소실 빠름
)

cyp2d6_labels <- c(
  "poor"   = "Poor metabolizer (8%)",
  "normal" = "Normal metabolizer (70%)",
  "rapid"  = "Rapid metabolizer (22%)"
)

df_plot <- results %>%
  filter(!is.na(t_50pct)) %>%   # v2: <1h 제외 규칙 삭제(D1/D4). censored/calc_failed는
                                  # t_50pct가 NA이므로 이 필터만으로 자동 제외되고,
                                  # reached 환자 분모는 Table 3/6에 별도 보고됨
  mutate(
    drug      = factor(drug, levels = drug_order),
    cyp2d6    = factor(cyp2d6, levels = c("poor", "normal", "rapid")),
    drug_label = factor(
      drug_labels[as.character(drug)],
      levels = unname(drug_labels[drug_order])
    )
  )

# ── 4. Median summary (console) ───────────────────────────────────────────────
cat("\n── Median t50% by drug & CYP2D6 (hours) ──\n")
df_plot %>%
  group_by(drug, cyp2d6) %>%
  summarise(
    median_h = round(median(t_50pct), 1),
    q25_h    = round(quantile(t_50pct, 0.25), 1),
    q75_h    = round(quantile(t_50pct, 0.75), 1),
    n        = n(),
    .groups  = "drop"
  ) %>%
  print(n = Inf)

# ── 5. Theme ──────────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor    = element_blank(),
    panel.grid.major.x  = element_blank(),       # x 격자 제거 — 깔끔
    panel.grid.major.y  = element_line(colour = "grey88", linewidth = 0.4),
    axis.title          = element_text(size = 10),
    axis.text           = element_text(size = 9),
    axis.text.x         = element_text(size = 9.5, lineheight = 1.25),
    legend.position     = "bottom",
    legend.title        = element_text(size = 9, face = "bold"),
    legend.text         = element_text(size = 9),
    legend.key.size     = unit(0.55, "cm"),
    legend.box.background = element_rect(colour = "grey80", linewidth = 0.3),
    plot.caption        = element_text(size = 8, colour = "grey40",
                                       hjust = 0, lineheight = 1.3),
    plot.margin         = margin(t = 0.3, r = 1.0, b = 1.2, l = 0.3, unit = "cm")
  )

# ── 6. Figure ─────────────────────────────────────────────────────────────────
p_fig3 <- ggplot(df_plot,
                 aes(x     = drug_label,
                     y     = t_50pct,
                     fill  = cyp2d6,
                     colour = cyp2d6)) +
  
  # 임상 참조선 (1일, 7일)
  geom_hline(yintercept = 24,  linetype = "dashed", colour = "grey60",
             linewidth = 0.4) +
  geom_hline(yintercept = 168, linetype = "dashed", colour = "grey60",
             linewidth = 0.4) +
  
  # Violin — 전체 분포
  geom_violin(
    position  = position_dodge(width = 0.85),
    alpha     = 0.30,
    linewidth = 0.35,
    trim      = TRUE,
    scale     = "width"    # 모든 violin 동일 최대 너비
  ) +
  
  # Boxplot overlay — IQR + 중앙값
  geom_boxplot(
    position      = position_dodge(width = 0.85),
    width         = 0.10,
    outlier.shape = NA,     # outlier 제거 (violin에서 이미 표현)
    fill          = "white",
    alpha         = 0.85,
    linewidth     = 0.40,
    colour        = "grey30"
  ) +
  
  # 참조선 라벨 (좌측 배치 — Fluoxetine violin 겹침 방지)
  annotate("text", x = 0.42, y = 24,
           label = "1 day", size = 2.6, colour = "grey50", hjust = 0) +
  annotate("text", x = 0.42, y = 168,
           label = "7 days", size = 2.6, colour = "grey50", hjust = 0) +
  
  # 색상 스케일
  scale_fill_manual(
    values = cyp2d6_colors,
    labels = cyp2d6_labels,
    name   = "CYP2D6 phenotype"
  ) +
  scale_colour_manual(
    values = cyp2d6_colors,
    labels = cyp2d6_labels,
    name   = "CYP2D6 phenotype"
  ) +
  
  # y축: log10 — 하한 1h, 상한 ~400d (v2: fluoxetine 관찰창 확장으로 범위 넓힘, D4)
  scale_y_log10(
    limits = c(1, 10000),
    breaks = c(1, 3, 6, 10, 24, 72, 168, 720, 2160, 4320, 9600),
    labels = c("1 h", "3 h", "6 h", "10 h", "1 d", "3 d",
               "7 d", "30 d", "90 d", "180 d", "400 d")
  ) +
  annotation_logticks(sides = "l", linewidth = 0.3, colour = "grey50") +
  
  # coord_cartesian: clip off for annotations, x limits only
  coord_cartesian(clip = "off", xlim = c(0.5, 5.5)) +
  
  # 범례 통합 (fill + colour 합치기)
  guides(fill   = guide_legend(title = "CYP2D6 phenotype",
                               override.aes = list(alpha = 0.7)),
         colour = guide_legend(title = "CYP2D6 phenotype")) +
  
  labs(
    x       = NULL,
    y       = expression(italic(t)[50]*"% C"[0]*"  (log scale)"),
    caption = "13 censored poor-metabolizer paroxetine patients (window ended before\nt50% C0 was reached) are excluded from this distribution; see Table 3/6 for denominators."
  ) +
  base_theme

# ── 7. Save ───────────────────────────────────────────────────────────────────
fig_path <- file.path(FIG_DIR, "figure3_mc_t50_violin_v2.png")
ggsave(fig_path, p_fig3, width = 12, height = 7, dpi = 300)
cat(sprintf("\n✓ Figure 3 saved: %s\n", fig_path))