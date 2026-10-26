# =============================================================================
# figure4_mc_dCdt_ridge_v2.R
# Figure 4 (v2): Monte Carlo — relative dC/dt @24h distribution (2-panel ridge)
#
# v2 changes: reads v2 MC output; column renamed dCdt_rel_24h -> reldiff_24h;
# Panel A x-range widened (0-0.015 /h) to fit fluoxetine's new tail under D1/D3.
#
# Panel A: Fluoxetine+NFX  (slow decline, x: 0–0.015 /h)
# Panel B: Sertraline, Paroxetine, Venlafaxine IR/XR  (fast decline)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)
  library(patchwork)
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

df_base <- results %>%
  filter(!is.na(reldiff_24h), reldiff_24h < 0) %>%   # v2: column renamed reldiff_24h
  mutate(abs_dCdt = abs(reldiff_24h),
         drug     = factor(drug))

# ── 2. Summary ────────────────────────────────────────────────────────────────
cat("\n── Median |dCdt_rel_24h| by drug ──\n")
df_base %>%
  group_by(drug) %>%
  summarise(median = round(median(abs_dCdt), 5),
            q25    = round(quantile(abs_dCdt, 0.25), 5),
            q75    = round(quantile(abs_dCdt, 0.75), 5),
            .groups = "drop") %>%
  print()

# ── 3. Shared theme ───────────────────────────────────────────────────────────
ridge_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.4),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9),
    axis.text.y        = element_text(size = 9.5),
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 11)
  )

x_label <- expression(
  group("|", frac(Delta*C[total], C[0]%.%Delta*t), "|") * "  (h"^{-1}*")"
)

# ── 4. Panel A: Fluoxetine ────────────────────────────────────────────────────
df_A <- df_base %>%
  filter(drug == "Fluoxetine") %>%
  mutate(drug_label = "Fluoxetine\n+NFX")

p_A <- ggplot(df_A, aes(x = abs_dCdt, y = drug_label,
                        fill = drug, colour = drug)) +
  geom_density_ridges(
    scale = 1.0, alpha = 0.65, linewidth = 0.45,
    rel_min_height = 0.005,
    quantile_lines = TRUE, quantiles = 0.5,
    vline_linetype = "dashed", vline_color = "grey20", vline_size = 0.6
  ) +
  scale_fill_manual(values   = drug_colors) +
  scale_colour_manual(values = drug_colors) +
  scale_x_continuous(
    limits = c(0, 0.015),   # v2: widened — reldiff_24h p99 ~0.0067, max ~0.0142
    breaks = c(0, 0.0025, 0.005, 0.0075, 0.010, 0.0125, 0.015),
    labels = c("0", "0.0025", "0.005", "0.0075", "0.010", "0.0125", "0.015"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  labs(title = "A.  Fluoxetine + NFX  (slow-decline group)",
       x = x_label, y = NULL) +
  ridge_theme

# ── 5. Panel B: Other 4 drugs ─────────────────────────────────────────────────
drug_order_B  <- c("Sertraline", "Paroxetine", "Venlafaxine_IR", "Venlafaxine_XR")
drug_labels_B <- c(
  "Sertraline"     = "Sertraline",
  "Paroxetine"     = "Paroxetine",
  "Venlafaxine_IR" = "Venlafaxine IR",
  "Venlafaxine_XR" = "Venlafaxine XR"
)

df_B <- df_base %>%
  filter(drug != "Fluoxetine") %>%
  mutate(drug_label = factor(drug_labels_B[as.character(drug)],
                             levels = unname(drug_labels_B[drug_order_B])))

x_upper_B <- quantile(df_B$abs_dCdt, 0.99, na.rm = TRUE)

p_B <- ggplot(df_B, aes(x = abs_dCdt, y = drug_label,
                        fill = drug, colour = drug)) +
  geom_density_ridges(
    scale = 1.8, alpha = 0.55, linewidth = 0.45,
    rel_min_height = 0.005,
    quantile_lines = TRUE, quantiles = 0.5,
    vline_linetype = "dashed", vline_color = "grey20", vline_size = 0.6
  ) +
  scale_fill_manual(values   = drug_colors) +
  scale_colour_manual(values = drug_colors) +
  scale_x_continuous(
    limits = c(0, x_upper_B),
    labels = scales::label_number(accuracy = 0.005),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  labs(title = "B.  SSRI / SNRI  (fast-decline group)",
       x = x_label, y = NULL) +
  ridge_theme

# ── 6. Assemble & save ────────────────────────────────────────────────────────
fig4 <- p_A + p_B +
  plot_layout(widths = c(1, 2.2))

fig_path <- file.path(FIG_DIR, "figure4_mc_dCdt_ridge_v2.png")
ggsave(fig_path, fig4, width = 13, height = 6, dpi = 300)
cat(sprintf("\n✓ Figure 4 saved: %s\n", fig_path))