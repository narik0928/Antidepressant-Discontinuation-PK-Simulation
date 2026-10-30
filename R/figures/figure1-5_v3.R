# =============================================================================
# figure1-5_v3.R
# Author: Nari Kim (code assist: Claude)
#
# Regenerates Figures 1-5 for TCP-26-013 manuscript v3 (round3_1029), using
# v3 parameters (data/pk_parameters/pk_parameters_v3.R) and v3 Monte Carlo
# output (outputs/mc_v3/mc_results_v3_all_n10000.csv).
#
# Adapted from the v2 scripts (R/figures/figure1_discontinuation_v2.R,
# figure2_venlafaxine_IR_XR_v2.R, figure3_mc_t50_violin_v2.R,
# figure4_mc_dCdt_ridge_v2.R, figure5_mc_scatter_v2.R); combined into one
# file per NARI_TASKS_1030.md Section 1.
#
# v3 changes from v2 relevant to these figures:
#   - Paroxetine dose is 30 mg/day (v2 and v1 used 20 mg/day) -- Figure 1C.
#   - CYP2D6 is modelled for paroxetine only; no rapid/ultrarapid category
#     exists in v3 (Methods 2.4.2). Figure 3 therefore splits ONLY the
#     paroxetine slot into poor/normal-metabolizer violins; the other four
#     drugs are plotted as single (non-split) violins. This differs from
#     v2's figure3, which dodge-split every drug by 3 phenotypes.
#   - Figure 4's x-axis is relabelled "mean relative decline rate over 0-24 h
#     (h^-1)" per NARI_TASKS_1030.md Section 1 (metric definition unchanged:
#     |C0 - C24| / (24 x C0); column `reldiff_24h` in the v3 MC output).
#   - Figure 2 legend values in the manuscript (IR C0=147.2 ng/mL,
#     t50%C0=15.2 h; XR C0=119.4, t50%C0=14.6 h) were computed by Dukhyun's
#     side from the v3 Table 1 typical parameter set (CHANGES_1029.md
#     Section 4). If this script's recomputed values differ from those by
#     more than 3%, STOP and contact Dukhyun rather than editing the
#     manuscript legend (NARI_TASKS_1030.md Section 1).
#
# Run order: source() this file as a whole, or run each numbered section
# (1-5) independently -- each section is self-contained after Section 0.
# =============================================================================

suppressPackageStartupMessages({
  library(rxode2)
  library(tidyverse)
  library(patchwork)
  library(ggridges)
  library(ggrepel)
})

# ── 0. Setup ──────────────────────────────────────────────────────────────────
BASE    <- "~/nari-research/pkpd-antidepressant-sim"
FIG_DIR <- file.path(BASE, "outputs", "figures")
MC_DIR  <- file.path(BASE, "outputs", "mc_v3")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "data", "pk_parameters", "pk_parameters_v3.R"))

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

run_sim <- function(model, params, inits, dose_F, tau_h, n_doses, obs_h, dt) {
  t_disc <- n_doses * tau_h
  ev <- et(amt = dose_F, ii = tau_h, addl = n_doses - 1, time = 0) %>%
    et(seq(0, t_disc + obs_h, by = dt))
  list(
    df      = rxSolve(model, params, ev, inits = inits) %>% as.data.frame(),
    t_disc  = t_disc,
    tau_h   = tau_h,
    n_doses = n_doses
  )
}

get_post <- function(res) {
  res$df %>%
    filter(time >= res$t_disc) %>%
    mutate(days = (time - res$t_disc) / 24)
}

# time (h, post-discontinuation) at which `conc` first falls to `frac` of C0
find_t_frac <- function(post_df, conc_col, frac) {
  C0 <- post_df[[conc_col]][1]
  target <- C0 * frac
  below <- which(post_df[[conc_col]] <= target)
  if (length(below) == 0) return(NA_real_)
  i <- below[1]
  if (i == 1) return(post_df$days[1] * 24)
  # linear interpolation between i-1 and i
  t0 <- post_df$days[i - 1] * 24; t1 <- post_df$days[i] * 24
  c0 <- post_df[[conc_col]][i - 1]; c1 <- post_df[[conc_col]][i]
  t0 + (c0 - target) / (c0 - c1) * (t1 - t0)
}

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.4),
    plot.title        = element_text(face = "bold", size = 11),
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 9),
    legend.background = element_rect(fill = alpha("white", 0.85),
                                     colour = "grey80", linewidth = 0.3),
    legend.key.width  = unit(1.5, "cm"),
    legend.text       = element_text(size = 9),
    legend.title      = element_blank()
  )

# =============================================================================
# 1. FIGURE 1 -- Post-discontinuation profiles, typical parameter set
#    Panels: A Sertraline | B Fluoxetine+NFX | C Paroxetine | D Venlafaxine XR+ODV
# =============================================================================

## 1A. Sertraline (50 mg QD x 14 doses; obs 30 d)
p_s <- pk_params$sertraline
res_sert <- run_sim(ode_1cmt,
  params = c(Ka = p_s$Ka, Vd = p_s$Vd, CL = p_s$CL),
  inits  = c(depot = 0, Cp = 0),
  dose_F = p_s$F * 50, tau_h = 24, n_doses = 14, obs_h = 30 * 24, dt = 1)
post_sert <- get_post(res_sert) %>% mutate(Cp_ng = Cp * 1000)
C0_sert <- post_sert$Cp_ng[1]

## 1B. Fluoxetine + Norfluoxetine (20 mg QD x 60 doses; obs 60 d)
p_f  <- pk_params$fluoxetine
nfx  <- pk_params$fluoxetine$norfluoxetine
res_flx <- run_sim(ode_parent_met,
  params = c(Ka = p_f$Ka, Vd_p = p_f$Vd, CL_p = p_f$CL,
             Fm = nfx$Fm, Vd_m = nfx$Vd_met, CL_m = nfx$CL_met),
  inits  = c(depot = 0, Cp = 0, Cm = 0),
  dose_F = p_f$F * 20, tau_h = 24, n_doses = 60, obs_h = 60 * 24, dt = 1)
post_flx <- get_post(res_flx)
post_flx_long <- post_flx %>%
  transmute(days,
    Fluoxetine      = Cp * 1000,
    Norfluoxetine   = Cm * 1000,
    `FLX+NFX total` = Ctotal * 1000) %>%
  pivot_longer(c(Fluoxetine, Norfluoxetine, `FLX+NFX total`),
               names_to = "component", values_to = "conc_ng") %>%
  mutate(component = factor(component,
         levels = c("FLX+NFX total", "Fluoxetine", "Norfluoxetine")))

## 1C. Paroxetine -- v3: 30 mg/day (D6), 14 doses, obs 14 d
p_p <- pk_params$paroxetine
res_par <- run_sim(ode_1cmt,
  params = c(Ka = p_p$Ka, Vd = p_p$Vd, CL = p_p$CL),
  inits  = c(depot = 0, Cp = 0),
  dose_F = p_p$F * 30, tau_h = 24, n_doses = 14, obs_h = 14 * 24, dt = 0.5)
post_par <- get_post(res_par) %>% mutate(Cp_ng = Cp * 1000)
C0_par <- post_par$Cp_ng[1]

## 1D. Venlafaxine XR + ODV (75 mg QD x 7 doses; obs 5 d)
p_v  <- pk_params$venlafaxine_XR
odv  <- pk_params$venlafaxine_XR$ODV
res_ven <- run_sim(ode_parent_met,
  params = c(Ka = p_v$Ka, Vd_p = p_v$Vd, CL_p = p_v$CL,
             Fm = odv$Fm, Vd_m = odv$Vd_odv, CL_m = odv$CL_odv),
  inits  = c(depot = 0, Cp = 0, Cm = 0),
  dose_F = p_v$F * 75, tau_h = 24, n_doses = 7, obs_h = 5 * 24, dt = 0.25)
post_ven <- get_post(res_ven)
post_ven_long <- post_ven %>%
  transmute(days,
    Venlafaxine     = Cp * 1000,
    ODV             = Cm * 1000,
    `VEN+ODV total` = Ctotal * 1000) %>%
  pivot_longer(c(Venlafaxine, ODV, `VEN+ODV total`),
               names_to = "component", values_to = "conc_ng") %>%
  mutate(component = factor(component,
         levels = c("VEN+ODV total", "Venlafaxine", "ODV")))

## Panels
p1A <- ggplot(post_sert, aes(days, Cp_ng)) +
  geom_line(colour = "#2166AC", linewidth = 0.8) +
  geom_hline(yintercept = C0_sert * c(0.75, 0.50), linetype = c("dashed", "dotted"),
             colour = "grey40", linewidth = 0.4) +
  labs(title = "A. Sertraline", x = "Days after discontinuation", y = "Concentration (ng/mL)") +
  base_theme

p1B <- ggplot(post_flx_long, aes(days, conc_ng, colour = component)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = post_flx$Ctotal[1] * 1000 * c(0.75, 0.50),
             linetype = c("dashed", "dotted"), colour = "grey40", linewidth = 0.4) +
  scale_colour_manual(values = c("FLX+NFX total" = "#2D6A4F",
                                 "Fluoxetine" = "#52B788", "Norfluoxetine" = "#95D5B2")) +
  labs(title = "B. Fluoxetine + Norfluoxetine", x = "Days after discontinuation",
       y = "Concentration (ng/mL)") +
  base_theme + theme(legend.position = c(0.8, 0.85))

p1C <- ggplot(post_par, aes(days, Cp_ng)) +
  geom_line(colour = "#D95F02", linewidth = 0.8) +
  geom_hline(yintercept = C0_par * c(0.75, 0.50), linetype = c("dashed", "dotted"),
             colour = "grey40", linewidth = 0.4) +
  labs(title = "C. Paroxetine", x = "Days after discontinuation", y = "Concentration (ng/mL)") +
  base_theme

p1D <- ggplot(post_ven_long, aes(days, conc_ng, colour = component)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = post_ven$Ctotal[1] * 1000 * c(0.75, 0.50),
             linetype = c("dashed", "dotted"), colour = "grey40", linewidth = 0.4) +
  scale_colour_manual(values = c("VEN+ODV total" = "#264653",
                                 "Venlafaxine" = "#5B9BD5", "ODV" = "#A9C9E8")) +
  labs(title = "D. Venlafaxine XR + ODV", x = "Days after discontinuation",
       y = "Concentration (ng/mL)") +
  base_theme + theme(legend.position = c(0.8, 0.85))

fig1 <- (p1A + p1B) / (p1C + p1D)
ggsave(file.path(FIG_DIR, "figure1_discontinuation_v3.jpg"), fig1, device = "jpeg",
       width = 12, height = 9, units = "in", dpi = 600, quality = 100)
cat(sprintf("\n✓ Figure 1 saved. Deterministic C0: sertraline %.1f ng/mL, paroxetine %.1f ng/mL\n",
            C0_sert, C0_par))

# =============================================================================
# 2. FIGURE 2 -- Venlafaxine IR vs XR, matched total daily dose (75 mg/day)
# =============================================================================

p_ir  <- pk_params$venlafaxine_IR
odv_ir <- pk_params$venlafaxine_IR$ODV
res_ir <- run_sim(ode_parent_met,
  params = c(Ka = p_ir$Ka, Vd_p = p_ir$Vd, CL_p = p_ir$CL,
             Fm = odv_ir$Fm, Vd_m = odv_ir$Vd_odv, CL_m = odv_ir$CL_odv),
  inits  = c(depot = 0, Cp = 0, Cm = 0),
  dose_F = p_ir$F * 25, tau_h = 8, n_doses = 21, obs_h = 7 * 24, dt = 0.1)
post_ir <- get_post(res_ir)
C0_ir   <- post_ir$Ctotal[1] * 1000
t50_ir  <- find_t_frac(post_ir %>% mutate(Ctotal_ng = Ctotal * 1000), "Ctotal_ng", 0.50)

p_xr  <- pk_params$venlafaxine_XR
odv_xr <- pk_params$venlafaxine_XR$ODV
res_xr <- run_sim(ode_parent_met,
  params = c(Ka = p_xr$Ka, Vd_p = p_xr$Vd, CL_p = p_xr$CL,
             Fm = odv_xr$Fm, Vd_m = odv_xr$Vd_odv, CL_m = odv_xr$CL_odv),
  inits  = c(depot = 0, Cp = 0, Cm = 0),
  dose_F = p_xr$F * 75, tau_h = 24, n_doses = 7, obs_h = 7 * 24, dt = 0.1)
post_xr <- get_post(res_xr)
C0_xr  <- post_xr$Ctotal[1] * 1000
t50_xr <- find_t_frac(post_xr %>% mutate(Ctotal_ng = Ctotal * 1000), "Ctotal_ng", 0.50)

cat(sprintf("\n— Figure 2 deterministic values (compare with manuscript legend: IR 147.2 ng/mL/15.2 h; XR 119.4 ng/mL/14.6 h) —\n"))
cat(sprintf("IR: C0 = %.1f ng/mL, t50%%C0 = %.1f h\n", C0_ir, t50_ir))
cat(sprintf("XR: C0 = %.1f ng/mL, t50%%C0 = %.1f h\n", C0_xr, t50_xr))
pct_diff_ir <- abs(C0_ir - 147.2) / 147.2 * 100
pct_diff_xr <- abs(C0_xr - 119.4) / 119.4 * 100
if (pct_diff_ir > 3 || pct_diff_xr > 3) {
  warning("Figure 2 recomputed C0 differs from the manuscript legend by more than 3% -- STOP, do not edit the manuscript, contact Dukhyun (NARI_TASKS_1030.md Section 1).")
}

post_ir_long <- post_ir %>%
  transmute(days, Venlafaxine = Cp * 1000, ODV = Cm * 1000, `VEN+ODV total` = Ctotal * 1000) %>%
  pivot_longer(c(Venlafaxine, ODV, `VEN+ODV total`), names_to = "component", values_to = "conc_ng") %>%
  mutate(component = factor(component, levels = c("VEN+ODV total", "Venlafaxine", "ODV")))
post_xr_long <- post_xr %>%
  transmute(days, Venlafaxine = Cp * 1000, ODV = Cm * 1000, `VEN+ODV total` = Ctotal * 1000) %>%
  pivot_longer(c(Venlafaxine, ODV, `VEN+ODV total`), names_to = "component", values_to = "conc_ng") %>%
  mutate(component = factor(component, levels = c("VEN+ODV total", "Venlafaxine", "ODV")))

y_max <- max(C0_ir, C0_xr) * 1.15
comp_colors <- c("VEN+ODV total" = "#264653", "Venlafaxine" = "#5B9BD5", "ODV" = "#A9C9E8")

p2A <- ggplot(post_ir_long, aes(days, conc_ng, colour = component)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = C0_ir * c(0.75, 0.50), linetype = c("dashed", "dotted"),
             colour = "grey40", linewidth = 0.4) +
  annotate("text", x = 6.8, y = y_max * 0.95, hjust = 1, size = 3.2,
           label = sprintf("C0 = %.1f ng/mL\nt50%%C0 = %.1f h", C0_ir, t50_ir)) +
  scale_colour_manual(values = comp_colors) +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(title = "A. Venlafaxine IR (25 mg q8h)", x = "Days after discontinuation",
       y = "Concentration (ng/mL)") +
  base_theme + theme(legend.position = "bottom")

p2B <- ggplot(post_xr_long, aes(days, conc_ng, colour = component)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = C0_xr * c(0.75, 0.50), linetype = c("dashed", "dotted"),
             colour = "grey40", linewidth = 0.4) +
  annotate("text", x = 6.8, y = y_max * 0.95, hjust = 1, size = 3.2,
           label = sprintf("C0 = %.1f ng/mL\nt50%%C0 = %.1f h", C0_xr, t50_xr)) +
  scale_colour_manual(values = comp_colors) +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(title = "B. Venlafaxine XR (75 mg QD)", x = "Days after discontinuation",
       y = "Concentration (ng/mL)") +
  base_theme + theme(legend.position = "bottom")

fig2 <- p2A + p2B + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "figure2_venlafaxine_IR_XR_v3.jpg"), fig2, device = "jpeg",
       width = 12, height = 6, units = "in", dpi = 600, quality = 100)

# =============================================================================
# 3. FIGURE 3 -- Monte Carlo t50% C0 distribution by drug; paroxetine split
#    by CYP2D6 phenotype only (v3: no rapid/ultrarapid category exists)
# =============================================================================

N <- 10000
results <- read_csv(file.path(MC_DIR, sprintf("mc_results_v3_all_n%d.csv", N)),
                     show_col_types = FALSE)

status_report <- results %>%
  group_by(drug, status_50pct) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(drug) %>%
  mutate(pct = round(n / sum(n) * 100, 2)) %>%
  ungroup()
cat("\n── t50% status report (reached / censored / calc_failed) ──\n")
print(status_report, n = Inf)

drug_order <- c("Venlafaxine_XR", "Venlafaxine_IR", "Paroxetine", "Sertraline", "Fluoxetine")
drug_labels <- c(
  "Venlafaxine_XR" = "Venlafaxine\nXR", "Venlafaxine_IR" = "Venlafaxine\nIR",
  "Paroxetine" = "Paroxetine", "Sertraline" = "Sertraline",
  "Fluoxetine" = "Fluoxetine\n+NFX")

df3 <- results %>%
  filter(!is.na(t_50pct)) %>%
  mutate(
    drug = factor(drug, levels = drug_order),
    drug_label = factor(drug_labels[as.character(drug)], levels = unname(drug_labels[drug_order])),
    # split field: paroxetine -> cyp2d6 (poor/normal); every other drug -> "n/a"
    # (single violin, not phenotype-dodged, since only paroxetine carries a
    # modelled CYP2D6 effect in v3 -- Methods 2.4.2)
    split = if_else(drug == "Paroxetine", cyp2d6, "n/a")
  )

cyp2d6_colors <- c("poor" = "#C0392B", "normal" = "#2980B9", "n/a" = "#6C6C6C")
cyp2d6_labels <- c("poor" = "Paroxetine, poor metabolizer (8%)",
                   "normal" = "Paroxetine, normal metabolizer (92%)",
                   "n/a" = "No CYP2D6 effect modelled")

p_fig3 <- ggplot(df3, aes(x = drug_label, y = t_50pct, fill = split, colour = split)) +
  geom_hline(yintercept = 24, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = 168, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_violin(position = position_dodge(width = 0.85), alpha = 0.30,
              linewidth = 0.35, trim = TRUE, scale = "width") +
  geom_boxplot(position = position_dodge(width = 0.85), width = 0.10,
               outlier.shape = NA, fill = "white", alpha = 0.85,
               linewidth = 0.40, colour = "grey30") +
  annotate("text", x = 0.42, y = 24, label = "1 day", size = 2.6, colour = "grey50", hjust = 0) +
  annotate("text", x = 0.42, y = 168, label = "7 days", size = 2.6, colour = "grey50", hjust = 0) +
  scale_fill_manual(values = cyp2d6_colors, name = NULL,
                     breaks = c("poor", "normal"),
                     labels = cyp2d6_labels[c("poor", "normal")]) +
  scale_colour_manual(values = cyp2d6_colors, name = NULL,
                       breaks = c("poor", "normal"),
                       labels = cyp2d6_labels[c("poor", "normal")]) +
  scale_y_log10(limits = c(1, 10000),
                breaks = c(1, 3, 6, 10, 24, 72, 168, 720, 2160, 4320, 9600),
                labels = c("1 h", "3 h", "6 h", "10 h", "1 d", "3 d", "7 d", "30 d", "90 d", "180 d", "400 d")) +
  annotation_logticks(sides = "l", linewidth = 0.3, colour = "grey50") +
  coord_cartesian(clip = "off", xlim = c(0.5, 5.5)) +
  labs(x = NULL, y = expression(italic(t)[50] * "% C"[0] * "  (log scale)")) +
  base_theme + theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "figure3_mc_t50_violin_v3.jpg"), p_fig3, device = "jpeg",
       width = 12, height = 7, units = "in", dpi = 600, quality = 100)
cat("\n✓ Figure 3 saved.\n")

# =============================================================================
# 4. FIGURE 4 -- Mean relative decline rate over 0-24 h after discontinuation
#    x-axis label updated per NARI_TASKS_1030.md Section 1
# =============================================================================

drug_colors <- c("Fluoxetine" = "#2D6A4F", "Sertraline" = "#2166AC",
                 "Paroxetine" = "#D95F02", "Venlafaxine_IR" = "#7B2D8B",
                 "Venlafaxine_XR" = "#264653")

df4_base <- results %>%
  filter(!is.na(reldiff_24h), reldiff_24h < 0) %>%
  mutate(abs_dCdt = abs(reldiff_24h), drug = factor(drug))

x_label4 <- "Mean relative decline rate over 0–24 h (h⁻¹)"

df4_A <- df4_base %>% filter(drug == "Fluoxetine") %>% mutate(drug_label = "Fluoxetine\n+NFX")
p4A <- ggplot(df4_A, aes(abs_dCdt, drug_label, fill = drug, colour = drug)) +
  geom_density_ridges(scale = 1.0, alpha = 0.65, linewidth = 0.45, rel_min_height = 0.005,
                       quantile_lines = TRUE, quantiles = 0.5,
                       vline_linetype = "dashed", vline_color = "grey20", vline_size = 0.6) +
  scale_fill_manual(values = drug_colors) + scale_colour_manual(values = drug_colors) +
  scale_x_continuous(limits = c(0, 0.005),
                      expand = expansion(mult = c(0, 0.02))) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  labs(title = "A. Fluoxetine + NFX (slow-decline group)", x = x_label4, y = NULL) +
  base_theme + theme(legend.position = "none", panel.grid.major.y = element_blank())

drug_order_B  <- c("Sertraline", "Paroxetine", "Venlafaxine_IR", "Venlafaxine_XR")
drug_labels_B <- c("Sertraline" = "Sertraline", "Paroxetine" = "Paroxetine",
                    "Venlafaxine_IR" = "Venlafaxine IR", "Venlafaxine_XR" = "Venlafaxine XR")
df4_B <- df4_base %>% filter(drug != "Fluoxetine") %>%
  mutate(drug_label = factor(drug_labels_B[as.character(drug)], levels = unname(drug_labels_B[drug_order_B])))
x_upper_B <- quantile(df4_B$abs_dCdt, 0.99, na.rm = TRUE)   # clip uppermost 1%, per Figure 4 legend

p4B <- ggplot(df4_B, aes(abs_dCdt, drug_label, fill = drug, colour = drug)) +
  geom_density_ridges(scale = 1.8, alpha = 0.55, linewidth = 0.45, rel_min_height = 0.005,
                       quantile_lines = TRUE, quantiles = 0.5,
                       vline_linetype = "dashed", vline_color = "grey20", vline_size = 0.6) +
  scale_fill_manual(values = drug_colors) + scale_colour_manual(values = drug_colors) +
  scale_x_continuous(limits = c(0, x_upper_B), labels = scales::label_number(accuracy = 0.005),
                      expand = expansion(mult = c(0, 0.02))) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  labs(title = "B. Sertraline, paroxetine, venlafaxine IR/XR (fast-decline group)",
       x = x_label4, y = NULL) +
  base_theme + theme(legend.position = "none", panel.grid.major.y = element_blank())

fig4 <- p4A + p4B + plot_layout(widths = c(1, 2.2))
ggsave(file.path(FIG_DIR, "figure4_mc_dCdt_ridge_v3.jpg"), fig4, device = "jpeg",
       width = 13, height = 6, units = "in", dpi = 600, quality = 100)
cat("\n✓ Figure 4 saved.\n")

# =============================================================================
# 5. FIGURE 5 -- Absolute AUC deficit (7 d) vs t50% C0, drug-level stress map
# =============================================================================

drug_labels_clean <- c("Fluoxetine" = "Fluoxetine +NFX", "Sertraline" = "Sertraline",
                        "Paroxetine" = "Paroxetine", "Venlafaxine_IR" = "Venlafaxine IR",
                        "Venlafaxine_XR" = "Venlafaxine XR")

df5 <- results %>%
  filter(!is.na(t_50pct), !is.na(AUC_deficit_7d), AUC_deficit_7d > 0) %>%
  mutate(log_t50 = log10(t_50pct), drug = factor(drug, levels = names(drug_colors)),
         drug_lab = drug_labels_clean[as.character(drug)])

median_df5 <- df5 %>%
  group_by(drug, drug_lab) %>%
  summarise(med_log_t50 = median(log_t50), med_AUC = median(AUC_deficit_7d), .groups = "drop")

p_fig5 <- ggplot(df5, aes(log_t50, AUC_deficit_7d, colour = drug, fill = drug)) +
  stat_ellipse(level = 0.95, type = "norm", geom = "polygon", alpha = 0.12, linewidth = 0.4) +
  stat_ellipse(level = 0.50, type = "norm", geom = "polygon", alpha = 0.20, linewidth = 0.0) +
  geom_point(data = median_df5, aes(med_log_t50, med_AUC), size = 4.5, shape = 21,
             stroke = 0.8, colour = "grey20", fill = drug_colors[as.character(median_df5$drug)],
             inherit.aes = FALSE) +
  geom_label_repel(data = median_df5, aes(med_log_t50, med_AUC, label = drug_lab),
                    size = 3.2, fontface = "bold", colour = drug_colors[as.character(median_df5$drug)],
                    fill = alpha("white", 0.85), label.size = 0.25, box.padding = 0.5,
                    min.segment.length = 0.2, inherit.aes = FALSE, seed = 2025) +
  scale_colour_manual(values = drug_colors) + scale_fill_manual(values = drug_colors) +
  scale_x_continuous(breaks = log10(c(10, 24, 72, 168, 720, 2160, 4320, 9600)),
                      labels = c("10 h", "1 d", "3 d", "7 d", "30 d", "90 d", "180 d", "400 d"),
                      expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  labs(x = expression(italic(t)[50] * "% C"[0] * "  (log scale)"),
       y = "AUC deficit, 7 days (mg/L·h)") +
  base_theme + theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "figure5_mc_scatter_v3.jpg"), p_fig5, device = "jpeg",
       width = 11, height = 7, units = "in", dpi = 600, quality = 100)
cat("\n✓ Figure 5 saved. All v3 figures written to outputs/figures/*_v3.jpg\n")
