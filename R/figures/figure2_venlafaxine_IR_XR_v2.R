# =============================================================================
# figure2_venlafaxine_IR_XR_v2.R  (v2: D3 F=1 source; D1 C0 definition unchanged — already correct)
# Figure 2: Venlafaxine IR vs XR — formulation-sensitive PK simulation
#           under matched discontinuation assumptions
#
# Panel layout (1×2, shared y-axis range):
#   A. Venlafaxine IR (TID) + ODV
#   B. Venlafaxine XR (QD)  + ODV
#
# Key message: identical PK parameters except Ka →
#   IR trough lower than XR → different C0 → different discontinuation stress
#
# Units: ng/mL | Output: outputs/figures/figure2_venlafaxine_IR_XR.png
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

# ── 1. ODE Model ──────────────────────────────────────────────────────────────
ode_parent_met <- rxode2({
  d/dt(depot) = -Ka  * depot
  d/dt(Cp)    =  Ka  * depot / Vd_p - (CL_p / Vd_p) * Cp
  d/dt(Cm)    =  Fm  * CL_p  * Cp   / Vd_m - (CL_m  / Vd_m) * Cm
  Ctotal      =  Cp + Cm
})

# ── 2. Helpers ────────────────────────────────────────────────────────────────
run_sim <- function(params, inits, dose_F, tau_h, n_doses, obs_h, dt = 0.25) {
  t_disc <- n_doses * tau_h
  ev <- et(amt = dose_F, ii = tau_h, addl = n_doses - 1, time = 0) %>%
    et(seq(0, t_disc + obs_h, by = dt))
  list(
    df      = rxSolve(ode_parent_met, params, ev, inits = inits) %>%
      as.data.frame(),
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

# t50%: 선형 보간으로 Ctotal이 C0의 50%에 도달하는 시간 (h)
find_t50_h <- function(days, conc_ng, C0) {
  target <- 0.50 * C0
  idx    <- which(conc_ng <= target)[1]
  if (is.na(idx) || idx == 1) return(NA_real_)
  t1 <- days[idx - 1] * 24; t2 <- days[idx] * 24
  c1 <- conc_ng[idx - 1];   c2 <- conc_ng[idx]
  t1 + (target - c1) * (t2 - t1) / (c2 - c1)
}

# ── 3. Simulations ────────────────────────────────────────────────────────────
p_ir  <- pk_params$venlafaxine_IR
p_xr  <- pk_params$venlafaxine_XR
odv   <- pk_params$venlafaxine_IR$ODV   # IR/XR ODV 파라미터 동일

shared_params <- c(
  Vd_p = p_ir$Vd,   CL_p = p_ir$CL,
  Fm   = odv$Fm,
  Vd_m = odv$Vd_odv, CL_m = odv$CL_odv
)
inits <- c(depot = 0, Cp = 0, Cm = 0)

# ── IR: 25mg TID → total daily dose = 75mg/day (XR 75mg QD와 일치)
# 제형 비교의 confound 제거: Ka만 다르고 total daily dose 동일
res_IR <- run_sim(
  params  = c(Ka = p_ir$Ka, shared_params),   # Ka = 1.50 /h
  inits   = inits,
  dose_F  = p_ir$F * 25,    # 0.45 × 25mg = 11.25mg (TID → 75mg/day total)
  tau_h   = 8,               # TID
  n_doses = 21,              # 7일 × 3회
  obs_h   = 5 * 24
)

# ── XR: QD, 7일 투여 ──
res_XR <- run_sim(
  params  = c(Ka = p_xr$Ka, shared_params),   # Ka = 0.25 /h
  inits   = inits,
  dose_F  = p_xr$F * 75,    # 0.45 × 75mg = 33.75mg
  tau_h   = 24,              # QD
  n_doses = 7,
  obs_h   = 5 * 24
)

# ── Post-discontinuation 데이터 ──
post_IR <- get_post(res_IR)
post_XR <- get_post(res_XR)

# ── C0 (중단 시점 농도) ──
C0_IR <- post_IR$Ctotal[1] * 1000
C0_XR <- post_XR$Ctotal[1] * 1000

# ── t50% C0 (hours) ──
t50_IR_h <- find_t50_h(post_IR$days, post_IR$Ctotal * 1000, C0_IR)
t50_XR_h <- find_t50_h(post_XR$days, post_XR$Ctotal * 1000, C0_XR)

cat(sprintf("IR: C0_total = %.1f ng/mL | t50%% C0 = %.1f h\n", C0_IR, t50_IR_h))
cat(sprintf("XR: C0_total = %.1f ng/mL | t50%% C0 = %.1f h\n", C0_XR, t50_XR_h))

# ── Long format ──
to_long <- function(post_df, formulation) {
  post_df %>%
    transmute(
      days,
      `VEN+ODV total`  = Ctotal * 1000,
      `Venlafaxine`    = Cp     * 1000,
      `ODV`            = Cm     * 1000
    ) %>%
    pivot_longer(
      cols      = c(`VEN+ODV total`, `Venlafaxine`, `ODV`),
      names_to  = "component",
      values_to = "conc_ng"
    ) %>%
    mutate(
      component   = factor(component,
                           levels = c("VEN+ODV total", "Venlafaxine", "ODV")),
      formulation = formulation
    )
}

long_IR <- to_long(post_IR, "IR")
long_XR <- to_long(post_XR, "XR")

# ── 4. Shared aesthetics ──────────────────────────────────────────────────────
# y축: 패널별 독립 (IR C0 ≠ XR C0 → 공유 y축이 Panel B 공간 낭비)
# C0 직접 annotate → 제형별 절대 농도 차이를 수치로 전달

# ── 5. Shared aesthetics ──────────────────────────────────────────────────────
ven_colors <- c(
  "VEN+ODV total" = "#264653",
  "Venlafaxine"   = "#7B2D8B",
  "ODV"           = "#F4A261"
)
ven_linetypes <- c(
  "VEN+ODV total" = "solid",
  "Venlafaxine"   = "solid",
  "ODV"           = "dashed"
)

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor     = element_blank(),
    panel.grid.major     = element_line(colour = "grey92", linewidth = 0.4),
    plot.title           = element_text(face = "bold", size = 11),
    axis.title           = element_text(size = 10),
    axis.text            = element_text(size = 9),
    legend.background    = element_rect(fill  = alpha("white", 0.85),
                                        colour = "grey80", linewidth = 0.3),
    legend.key.width     = unit(1.5, "cm"),
    legend.text          = element_text(size = 9),
    legend.title         = element_blank()
  )

obs_days <- 5
x_breaks <- seq(0, obs_days, by = 1)

# ── 6. Build Panels ───────────────────────────────────────────────────────────

make_panel <- function(long_df, C0, t50_h, formulation_label, panel_letter,
                       show_y = TRUE, show_legend = FALSE) {
  
  # 패널별 독립 y축 범위
  y_max <- C0 * 1.20
  
  # t50% annotate 문자열
  t50_str <- if (!is.na(t50_h)) {
    sprintf("t\u2085\u2080%%\u00a0C\u2080\u00a0=\u00a0%.1f\u00a0h", t50_h)
  } else {
    "t\u2085\u2080%\u00a0C\u2080\u00a0>\u00a0120\u00a0h"
  }
  
  # C0 annotate 문자열
  c0_str <- sprintf("C\u2080\u00a0=\u00a0%.1f\u00a0ng/mL", C0)
  
  p <- ggplot(long_df, aes(x = days, y = conc_ng,
                           colour = component, linetype = component)) +
    # 임계선
    geom_hline(yintercept = 0.75 * C0,
               linetype = "dashed",  colour = "grey55", linewidth = 0.45) +
    geom_hline(yintercept = 0.50 * C0,
               linetype = "dotted", colour = "grey55", linewidth = 0.45) +
    # 농도 곡선 (VEN+ODV total 더 굵게)
    geom_line(aes(linewidth = component)) +
    # 임계선 라벨 (우측)
    annotate("text", x = obs_days * 0.97, y = 0.75 * C0 * 1.09,
             label = "75% C0", size = 2.8, colour = "grey40", hjust = 1) +
    annotate("text", x = obs_days * 0.97, y = 0.50 * C0 * 1.09,
             label = "50% C0", size = 2.8, colour = "grey40", hjust = 1) +
    # t50% 수치 (좌하단 — 곡선과 겹치지 않도록 y 낮춤)
    annotate("text", x = 0.20, y = y_max * 0.03,
             label = t50_str, size = 3.2, colour = "#264653",
             hjust = 0, fontface = "bold") +
    # C0 수치 (우상단 곡선 시작점 근처)
    annotate("text", x = 0.12, y = C0 * 1.08,
             label = c0_str, size = 3.0, colour = "#264653",
             hjust = 0, fontface = "italic") +
    scale_colour_manual(values   = ven_colors) +
    scale_linetype_manual(values = ven_linetypes) +
    scale_linewidth_manual(
      values = c("VEN+ODV total" = 1.4, "Venlafaxine" = 1.0, "ODV" = 1.0),
      guide  = "none"
    ) +
    scale_x_continuous(breaks = x_breaks, limits = c(0, obs_days)) +
    scale_y_continuous(limits = c(0, y_max),
                       expand = expansion(mult = c(0, 0))) +
    labs(
      title = sprintf("%s.  Venlafaxine %s  (SNRI)  + ODV",
                      panel_letter, formulation_label),
      x     = "Days after discontinuation",
      y     = if (show_y) "Plasma concentration (ng/mL)" else NULL
    ) +
    base_theme +
    theme(legend.position = if (show_legend) c(0.97, 0.97) else "none",
          legend.justification = c(1, 1))
  
  p
}

p_A <- make_panel(long_IR, C0_IR, t50_IR_h,
                  formulation_label = "IR",
                  panel_letter      = "A",
                  show_y            = TRUE,
                  show_legend       = FALSE)

p_B <- make_panel(long_XR, C0_XR, t50_XR_h,
                  formulation_label = "XR",
                  panel_letter      = "B",
                  show_y            = FALSE,
                  show_legend       = TRUE)

# ── 7. Assemble & Save ────────────────────────────────────────────────────────
fig2 <- p_A + p_B

fig_path <- file.path(FIG_DIR, "figure2_venlafaxine_IR_XR_v2.png")
ggsave(fig_path, fig2, width = 12, height = 5.5, dpi = 300)
cat(sprintf("\n✓ Figure 2 saved: %s\n", fig_path))