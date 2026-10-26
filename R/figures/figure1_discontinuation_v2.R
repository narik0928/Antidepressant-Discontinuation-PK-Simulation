# =============================================================================
# figure1_discontinuation_v2.R  (v2: D3 F=1 source; D1 C0 definition unchanged — already correct)
# Figure 1: Post-discontinuation plasma concentration profiles
#           4 drugs — 2×2 panel | Abrupt discontinuation
#           Population mean parameters | Operational PK thresholds: 75%, 50% Css
#
# Panel layout:
#   A. Sertraline (SSRI, IR)          B. Fluoxetine (SSRI, IR) + Norfluoxetine
#   C. Paroxetine (SSRI, IR)          D. Venlafaxine XR (SNRI) + ODV
#
# Units: ng/mL (simulation in mg/L × 1000)
# Output: outputs/figures/figure1_discontinuation.png
# =============================================================================

suppressPackageStartupMessages({
  library(rxode2)
  library(tidyverse)
  library(patchwork)
})

# ── 0. Setup ──────────────────────────────────────────────────────────────────
BASE    <- "~/nari-research/pkpd-antidepressant-sim"
FIG_DIR <- file.path(BASE, "outputs", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "data", "pk_parameters", "pk_parameters_v2.R"))

# ── 1. ODE Models ─────────────────────────────────────────────────────────────
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

# ── 2. Simulation Helper ──────────────────────────────────────────────────────
# Returns full simulation df + t_disc metadata
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

# Css: mean over last dosing interval (multi-dose SS)
calc_css <- function(res, conc_col) {
  df <- res$df
  t0 <- (res$n_doses - 1) * res$tau_h
  t1 <- res$t_disc
  mean(df[[conc_col]][df$time >= t0 & df$time <= t1], na.rm = TRUE)
}

# Post-discontinuation subset: re-zero time to days
get_post <- function(res) {
  res$df %>%
    filter(time >= res$t_disc) %>%
    mutate(days = (time - res$t_disc) / 24)
}

# ── 3. Shared Theme ───────────────────────────────────────────────────────────
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

# Threshold annotation helper
# Places label at x = 2% from right edge, slightly above threshold line
add_thresh_label <- function(p, obs_days, thresh_val, label_str, y_nudge = 0.08) {
  p + annotate(
    "text",
    x      = obs_days * 0.97,
    y      = thresh_val * (1 + y_nudge),
    label  = label_str,
    size   = 2.8,
    colour = "grey40",
    hjust  = 1
  )
}

# ── 4. Simulations ────────────────────────────────────────────────────────────

# ── 4A. Sertraline ────────────────────────────────────────────────────────────
p_s <- pk_params$sertraline
res_sert <- run_sim(
  model   = ode_1cmt,
  params  = c(Ka = p_s$Ka, Vd = p_s$Vd, CL = p_s$CL),
  inits   = c(depot = 0, Cp = 0),
  dose_F  = p_s$F * 50,     # 1.0 × 50mg
  tau_h   = 24,
  n_doses = 14,
  obs_h   = 30 * 24,
  dt      = 1
)
# C0 = 중단 시점 농도 — 임계선 기준 (mean Css 아님)
post_sert <- get_post(res_sert) %>% mutate(Cp_ng = Cp * 1000)
C0_sert   <- post_sert$Cp_ng[1]

# ── 4B. Fluoxetine + Norfluoxetine ────────────────────────────────────────────
p_f <- pk_params$fluoxetine
nfx <- pk_params$fluoxetine$norfluoxetine
res_flx <- run_sim(
  model   = ode_parent_met,
  params  = c(Ka   = p_f$Ka,
              Vd_p = p_f$Vd,   CL_p = p_f$CL,
              Fm   = nfx$Fm,
              Vd_m = nfx$Vd_met, CL_m = nfx$CL_met),
  inits   = c(depot = 0, Cp = 0, Cm = 0),
  dose_F  = p_f$F * 20,     # 0.70 × 20mg
  tau_h   = 24,
  n_doses = 60,
  obs_h   = 60 * 24,
  dt      = 1                # coarser dt for 120-day simulation
)
Css_flx_total <- calc_css(res_flx, "Ctotal") * 1000
post_flx <- get_post(res_flx)
C0_flx_total  <- (post_flx$Ctotal[1]) * 1000

# Long format for multi-line plot
post_flx_long <- post_flx %>%
  transmute(
    days,
    Fluoxetine    = Cp     * 1000,
    Norfluoxetine = Cm     * 1000,
    `FLX+NFX total` = Ctotal * 1000
  ) %>%
  pivot_longer(
    cols      = c(Fluoxetine, Norfluoxetine, `FLX+NFX total`),
    names_to  = "component",
    values_to = "conc_ng"
  ) %>%
  mutate(component = factor(component,
                            levels = c("FLX+NFX total", "Fluoxetine", "Norfluoxetine")))

# ── 4C. Paroxetine ────────────────────────────────────────────────────────────
p_p <- pk_params$paroxetine
res_par <- run_sim(
  model   = ode_1cmt,
  params  = c(Ka = p_p$Ka, Vd = p_p$Vd, CL = p_p$CL),
  inits   = c(depot = 0, Cp = 0),
  dose_F  = p_p$F * 20,     # 0.50 × 20mg
  tau_h   = 24,
  n_doses = 14,
  obs_h   = 14 * 24,
  dt      = 0.5
)
Css_par  <- calc_css(res_par, "Cp") * 1000
post_par <- get_post(res_par) %>% mutate(Cp_ng = Cp * 1000)
C0_par   <- post_par$Cp_ng[1]

# ── 4D. Venlafaxine XR + ODV ──────────────────────────────────────────────────
p_v  <- pk_params$venlafaxine_XR
odv  <- pk_params$venlafaxine_XR$ODV
res_ven <- run_sim(
  model   = ode_parent_met,
  params  = c(Ka   = p_v$Ka,
              Vd_p = p_v$Vd,   CL_p = p_v$CL,
              Fm   = odv$Fm,
              Vd_m = odv$Vd_odv, CL_m = odv$CL_odv),
  inits   = c(depot = 0, Cp = 0, Cm = 0),
  dose_F  = p_v$F * 75,     # 0.45 × 75mg
  tau_h   = 24,
  n_doses = 7,
  obs_h   = 5 * 24,
  dt      = 0.25
)
Css_ven_total <- calc_css(res_ven, "Ctotal") * 1000
post_ven <- get_post(res_ven)
C0_ven_total  <- (post_ven$Ctotal[1]) * 1000

post_ven_long <- post_ven %>%
  transmute(
    days,
    `Venlafaxine` = Cp     * 1000,
    `ODV`         = Cm     * 1000,
    `VEN+ODV total` = Ctotal * 1000
  ) %>%
  pivot_longer(
    cols      = c(`Venlafaxine`, `ODV`, `VEN+ODV total`),
    names_to  = "component",
    values_to = "conc_ng"
  ) %>%
  mutate(component = factor(component,
                            levels = c("VEN+ODV total", "Venlafaxine", "ODV")))

# ── 5. Build Panels ───────────────────────────────────────────────────────────

## ── Panel A: Sertraline ───────────────────────────────────────────────────────
obs_sert <- 30
p_A <- ggplot(post_sert, aes(x = days, y = Cp_ng)) +
  geom_hline(yintercept = 0.75 * C0_sert,
             linetype = "dashed",  colour = "grey55", linewidth = 0.45) +
  geom_hline(yintercept = 0.50 * C0_sert,
             linetype = "dotted", colour = "grey55", linewidth = 0.45) +
  geom_line(colour = "#2166AC", linewidth = 1.0) +
  annotate("text", x = obs_sert * 0.97, y = 0.75 * C0_sert * 1.10,
           label = "75% C0", size = 2.8, colour = "grey40", hjust = 1) +
  annotate("text", x = obs_sert * 0.97, y = 0.50 * C0_sert * 1.10,
           label = "50% C0", size = 2.8, colour = "grey40", hjust = 1) +
  scale_x_continuous(breaks = seq(0, obs_sert, by = 10),
                     limits = c(0, obs_sert)) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(title = "A.  Sertraline  (SSRI · IR)",
       x = NULL,
       y = "Plasma concentration (ng/mL)") +
  base_theme +
  theme(legend.position = "none")

## ── Panel B: Fluoxetine + Norfluoxetine ──────────────────────────────────────
obs_flx <- 60
flx_colors <- c(
  "FLX+NFX total" = "#2D6A4F",
  "Fluoxetine"    = "#E63946",
  "Norfluoxetine" = "#457B9D"
)
flx_lines <- c(
  "FLX+NFX total" = "solid",
  "Fluoxetine"    = "solid",
  "Norfluoxetine" = "dashed"
)

p_B <- ggplot(post_flx_long, aes(x = days, y = conc_ng,
                                 colour = component, linetype = component)) +
  geom_hline(yintercept = 0.75 * C0_flx_total,
             linetype = "dashed",  colour = "grey55", linewidth = 0.45) +
  geom_hline(yintercept = 0.50 * C0_flx_total,
             linetype = "dotted", colour = "grey55", linewidth = 0.45) +
  geom_line(linewidth = 0.9) +
  # 라벨을 중앙-우측에 배치 — 초기 곡선 겹침 방지
  annotate("text", x = 20, y = 0.75 * C0_flx_total * 1.10,
           label = "75% C0", size = 2.8, colour = "grey40", hjust = 0) +
  annotate("text", x = 20, y = 0.50 * C0_flx_total * 1.10,
           label = "50% C0", size = 2.8, colour = "grey40", hjust = 0) +
  scale_colour_manual(values = flx_colors) +
  scale_linetype_manual(values = flx_lines) +
  scale_x_continuous(breaks = seq(0, obs_flx, by = 15),
                     limits = c(0, obs_flx)) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(title = "B.  Fluoxetine  (SSRI · IR)  + Norfluoxetine",
       x = NULL,
       y = NULL) +
  base_theme +
  theme(legend.position  = c(0.97, 0.97),
        legend.justification = c(1, 1))

## ── Panel C: Paroxetine ───────────────────────────────────────────────────────
obs_par <- 14
p_C <- ggplot(post_par, aes(x = days, y = Cp_ng)) +
  geom_hline(yintercept = 0.75 * C0_par,
             linetype = "dashed",  colour = "grey55", linewidth = 0.45) +
  geom_hline(yintercept = 0.50 * C0_par,
             linetype = "dotted", colour = "grey55", linewidth = 0.45) +
  geom_line(colour = "#D95F02", linewidth = 1.0) +
  annotate("text", x = obs_par * 0.97, y = 0.75 * C0_par * 1.10,
           label = "75% C0", size = 2.8, colour = "grey40", hjust = 1) +
  annotate("text", x = obs_par * 0.97, y = 0.50 * C0_par * 1.10,
           label = "50% C0", size = 2.8, colour = "grey40", hjust = 1) +
  scale_x_continuous(breaks = seq(0, obs_par, by = 2),
                     limits = c(0, obs_par)) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(title = "C.  Paroxetine  (SSRI · IR)",
       x = "Days after discontinuation",
       y = "Plasma concentration (ng/mL)") +
  base_theme +
  theme(legend.position = "none")

## ── Panel D: Venlafaxine XR + ODV ─────────────────────────────────────────────
obs_ven <- 5
ven_colors <- c(
  "VEN+ODV total" = "#264653",
  "Venlafaxine"   = "#7B2D8B",
  "ODV"           = "#F4A261"
)
ven_lines <- c(
  "VEN+ODV total" = "solid",
  "Venlafaxine"   = "solid",
  "ODV"           = "dashed"
)

p_D <- ggplot(post_ven_long, aes(x = days, y = conc_ng,
                                 colour = component, linetype = component)) +
  geom_hline(yintercept = 0.75 * C0_ven_total,
             linetype = "dashed",  colour = "grey55", linewidth = 0.45) +
  geom_hline(yintercept = 0.50 * C0_ven_total,
             linetype = "dotted", colour = "grey55", linewidth = 0.45) +
  geom_line(linewidth = 1.1) +   # 굵기 증가 — VEN parent 선 가시성 확보
  annotate("text", x = obs_ven * 0.97, y = 0.75 * C0_ven_total * 1.10,
           label = "75% C0", size = 2.8, colour = "grey40", hjust = 1) +
  annotate("text", x = obs_ven * 0.97, y = 0.50 * C0_ven_total * 1.10,
           label = "50% C0", size = 2.8, colour = "grey40", hjust = 1) +
  scale_colour_manual(values = ven_colors) +
  scale_linetype_manual(values = ven_lines) +
  scale_x_continuous(breaks = seq(0, obs_ven, by = 1),
                     limits = c(0, obs_ven)) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(title = "D.  Venlafaxine XR  (SNRI)  + ODV",
       x = "Days after discontinuation",
       y = NULL) +
  base_theme +
  theme(legend.position  = c(0.97, 0.32),
        legend.justification = c(1, 1))

# ── 6. Assemble & Save ────────────────────────────────────────────────────────
fig1 <- (p_A + p_B) / (p_C + p_D)

fig_path <- file.path(FIG_DIR, "figure1_discontinuation_v2.png")
ggsave(fig_path, fig1, width = 12, height = 8, dpi = 300)
cat(sprintf("✓ Figure 1 saved: %s\n", fig_path))