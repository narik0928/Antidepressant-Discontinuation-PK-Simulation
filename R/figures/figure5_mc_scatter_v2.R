# =============================================================================
# figure5_mc_scatter_v2.R
# Figure 5 (v2): AUC deficit (7d, absolute mg/L*h) vs t50% C0 — drug-level PK stress map
#
# Geometry: 95% density ellipse per drug + median point + drug label
# NOTE: stat_ellipse on log10(t_50pct) to avoid distortion on log x-axis
# v2 changes: reads v2 MC output; no t_50pct >= 1 exclusion (NA = censored/
# calc_failed, dropped by !is.na alone); x-axis breaks extended to ~400 days
# for fluoxetine's widened observation window (D1/D3/D4). y-axis still shows
# the ABSOLUTE AUC deficit (matches Table 3); the normalised AUC deficit (%)
# used in Table 5/Section 3.7 is not plotted here — flagged for Nari/Dukhyun
# to decide whether a second panel or a swapped y-axis is wanted.
#
# Input : outputs/mc_v2/mc_results_v2_all_n10000.csv
# Output: outputs/figures/figure5_mc_scatter_v2.png
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

BASE    <- "~/nari-research/pkpd-antidepressant-sim"
FIG_DIR <- file.path(BASE, "outputs", "figures")
OUT_DIR <- file.path(BASE, "outputs", "mc_v2")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

N <- 10000

# ── 1. Load & prep ────────────────────────────────────────────────────────────
results <- read_csv(
  file.path(OUT_DIR, sprintf("mc_results_v2_all_n%d.csv", N)),
  show_col_types = FALSE
)

drug_colors <- c(
  "Fluoxetine"     = "#2D6A4F",
  "Sertraline"     = "#2166AC",
  "Paroxetine"     = "#D95F02",
  "Venlafaxine_IR" = "#7B2D8B",
  "Venlafaxine_XR" = "#264653"
)

drug_labels_clean <- c(
  "Fluoxetine"     = "Fluoxetine +NFX",
  "Sertraline"     = "Sertraline",
  "Paroxetine"     = "Paroxetine",
  "Venlafaxine_IR" = "Venlafaxine IR",
  "Venlafaxine_XR" = "Venlafaxine XR"
)

df_plot <- results %>%
  filter(
    !is.na(t_50pct),                      # v2: <1h 제외 규칙 삭제(D1/D4)
    !is.na(AUC_deficit_7d), AUC_deficit_7d > 0
  ) %>%
  mutate(
    log_t50  = log10(t_50pct),      # log 변환 — ellipse 계산 정확도
    drug     = factor(drug, levels = names(drug_colors)),
    drug_lab = drug_labels_clean[as.character(drug)]
  )

# ── 2. Median per drug ────────────────────────────────────────────────────────
median_df <- df_plot %>%
  group_by(drug, drug_lab) %>%
  summarise(
    med_log_t50 = median(log_t50),
    med_AUC     = median(AUC_deficit_7d),
    .groups     = "drop"
  )

cat("\n── Median values by drug ──\n")
median_df %>%
  mutate(med_t50_h = round(10^med_log_t50, 1),
         med_AUC   = round(med_AUC, 3)) %>%
  select(drug, med_t50_h, med_AUC) %>%
  print()

# ── 3. Theme ──────────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey90", linewidth = 0.4),
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 9),
    legend.position   = "none",
    plot.caption      = element_text(size = 8, colour = "grey40",
                                     hjust = 0, lineheight = 1.3),
    plot.margin       = margin(0.3, 0.5, 0.8, 0.3, "cm")
  )

# ── 4. Figure ─────────────────────────────────────────────────────────────────
p_fig5 <- ggplot(df_plot, aes(x = log_t50, y = AUC_deficit_7d,
                              colour = drug, fill = drug)) +
  
  # 95% density ellipse (log10 공간에서 계산 → 정확)
  stat_ellipse(
    level     = 0.95,
    type      = "norm",
    geom      = "polygon",
    alpha     = 0.12,
    linewidth = 0.4
  ) +
  stat_ellipse(
    level     = 0.50,
    type      = "norm",
    geom      = "polygon",
    alpha     = 0.20,
    linewidth = 0.0
  ) +
  
  # 중앙값 포인트
  geom_point(
    data     = median_df,
    mapping  = aes(x = med_log_t50, y = med_AUC),
    size     = 4.5,
    shape    = 21,
    stroke   = 0.8,
    colour   = "grey20",
    fill     = drug_colors[as.character(median_df$drug)],
    inherit.aes = FALSE
  ) +
  
  # 약물 라벨 (겹침 방지)
  geom_label_repel(
    data        = median_df,
    mapping     = aes(x = med_log_t50, y = med_AUC, label = drug_lab),
    size        = 3.2,
    fontface    = "bold",
    colour      = drug_colors[as.character(median_df$drug)],
    fill        = alpha("white", 0.85),
    label.size  = 0.25,
    box.padding = 0.5,
    min.segment.length = 0.2,
    inherit.aes = FALSE,
    seed        = 2025
  ) +
  
  scale_colour_manual(values = drug_colors) +
  scale_fill_manual(values   = drug_colors) +
  
  # x축: log10 공간 → 원래 단위 label
  scale_x_continuous(
    breaks = log10(c(5, 10, 24, 72, 168, 720, 2160, 4320, 9600)),
    labels = c("5 h", "10 h", "1 d", "3 d", "7 d", "30 d", "90 d", "180 d", "400 d"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.10))
  ) +
  
  labs(
    x       = expression(italic(t)[50]*"% C"[0]*"  (log scale)"),
    y       = "AUC deficit, 7 days  (mg/L·h)"
  ) +
  base_theme

# ── 5. Save ───────────────────────────────────────────────────────────────────
fig_path <- file.path(FIG_DIR, "figure5_mc_scatter_v2.png")
ggsave(fig_path, p_fig5, width = 11, height = 7, dpi = 300)
cat(sprintf("\n✓ Figure 5 saved: %s\n", fig_path))