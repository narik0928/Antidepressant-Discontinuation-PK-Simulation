# Antidepressant Discontinuation Risk — PK/PD Simulation Study

**Authors:** Nari Kim, Dukhyun Kim
**Institution:** Kingston University London, BSc Pharmacology
**Target Journal:** Translational and Clinical Pharmacology (TCP), manuscript TCP-26-013
**Status:** Resubmission (v2 re-analysis complete, 2026-10-26) — see `04_resubmission/`

---

## Research Overview

This repository contains the computational pharmacokinetic (PK) simulation
framework developed to quantify exposure-loss stress across selected SSRI
and SNRI antidepressants under abrupt discontinuation scenarios.

### Research Question
> Can pharmacokinetic properties — including half-life, active metabolite
> accumulation, formulation, and CYP2D6 phenotype — explain differential
> exposure-loss stress across selected SSRIs and SNRIs under abrupt
> discontinuation?

### v2 re-analysis (current, reported in the resubmitted manuscript)

During peer review the authors identified and corrected three errors not
raised by the reviewer, in addition to the reviewer-requested changes
(full detail in `04_resubmission/WORK_LOG.md` and
`04_resubmission/NUMBER_SOURCES.md`):

- **D1** — the discontinuation reference concentration (C0) was redefined
  as the instantaneous concentration at the moment of discontinuation
  (previously an interval-averaged concentration), matching the
  manuscript's stated definition.
- **D2** — PK parameter *values* were kept, but provenance labelling was
  corrected (author-estimated vs. literature-sourced vs. unconfirmed).
- **D3** — bioavailability (F) is now applied exactly once per regimen
  (F = 1.0 unified across all five regimens); it had previously been
  applied twice for four of the five regimens, inflating/deflating
  absolute concentrations by a per-drug constant factor.
- **D4** — the fluoxetine observation window was extended from 60 days to
  ~400 days, with explicit reached/censored/calc_failed status coding per
  simulated patient (the `< 1h` exclusion rule used in v1 was removed
  entirely).

### Key Findings (v2, corrected)

- Venlafaxine XR shows the fastest relative decline (median t50% C0 = 12.8
  [9.4–17.4] h, stress rank 1)
- Fluoxetine shows the slowest decline (median t50% C0 = 368.4
  [249.6–564.6] h) due to norfluoxetine buffering — a 28.7-fold difference
  vs. venlafaxine XR
- Simulated t50% C0 diverges from label-reported parent terminal half-life
  by 1.01–2.64-fold; the divergence is concentrated in the two
  active-metabolite agents (fluoxetine, venlafaxine), not the two
  parent-only agents (sertraline, paroxetine)
- Paroxetine shows the largest CYP2D6 effect (η² = 0.22, large), followed
  by fluoxetine (η² = 0.091), venlafaxine IR (η² = 0.063) and venlafaxine
  XR (η² = 0.055)
- A normalised AUC deficit metric (% of hypothetical constant exposure
  lost) separates drugs that the absolute AUC deficit could not
  statistically distinguish (fluoxetine vs. paroxetine)
- Censoring was negligible overall (13 of 50,000 simulated patients
  overall, all paroxetine poor-metabolizers who had not reached the t50%
  threshold within the observation window)

**v1 (original submission) numbers above have been superseded — do not
cite figures from the original `outputs/mc/`, `outputs/tables/`,
`outputs/qualification/` or `outputs/mc_sensitivity/` directories (non-`_v2`
paths) as current results. Those are kept only for audit/comparison
purposes (see `outputs/tables_v2/s1_v2_old_vs_new_paired.csv`).**

---

## Repository Structure

```
.
├── R/
│   ├── models/                     # Phase 3: Deterministic PK models (v1, unchanged)
│   │   ├── sertraline_model.R
│   │   ├── fluoxetine_model.R
│   │   ├── paroxetine_model.R
│   │   └── venlafaxine_model.R
│   │
│   ├── qualification/
│   │   ├── model_qualification.R   # v1 (superseded)
│   │   └── qualification_v2.R      # v2 — parent/metabolite Css reported separately
│   │
│   ├── simulation/
│   │   ├── mc_params.R             # Parameter sampling (n=10,000/drug, shared v1/v2 cohort)
│   │   ├── mc_run_single.R         # v1 (superseded)
│   │   └── mc_run_v2.R             # v2 — C0 redefinition (D1), F=1 (D3), status codes, extended fluoxetine window (D4)
│   │
│   ├── sensitivity/
│   │   ├── sensitivity_analysis.R  # v1 (superseded)
│   │   └── sensitivity_v2.R        # v2 — normalised vs. absolute AUC deficit split (A2 fix)
│   │
│   ├── analysis/
│   │   ├── stability_check.R       # n=1,000/5,000/10,000 convergence (unchanged)
│   │   ├── phase5_statistics.R     # v1 (superseded)
│   │   └── make_tables_v2.R        # v2 — Tables 3-6, S1, B1, Dunn tests from v2 output
│   │
│   └── figures/
│       ├── figure1_discontinuation.R / *_v2.R
│       ├── figure2_venlafaxine_IR_XR.R / *_v2.R
│       ├── figure3_mc_t50_violin.R / *_v2.R
│       ├── figure4_mc_dCdt_ridge.R / *_v2.R
│       └── figure5_mc_scatter.R / *_v2.R
│       (v1 scripts kept for reference; `_v2.R` versions generate the
│       figures reported in the resubmitted manuscript)
│
├── data/
│   └── pk_parameters/
│       ├── pk_parameters.R         # v1 (superseded)
│       └── pk_parameters_v2.R      # v2 — same parameter values as v1, F=1.0 unified, provenance labels corrected (D2/D3)
│
├── outputs/
│   ├── figures/                    # v1 + v2 PNGs (see naming: *_v2.png is current)
│   ├── qualification/ , qualification_v2/
│   ├── mc/ , mc_v2/
│   ├── mc_sensitivity/ , mc_sensitivity_v2/
│   └── tables/ , tables_v2/        # tables_v2/ is current; tables/ kept for audit
│
```

Manuscript drafts, the response letter, cover letter, and the full change
record (`WORK_LOG.md`, `NUMBER_SOURCES.md`, `RUN_ORDER.md`,
`_baseline_1026/` frozen pre-v2 docx backups) live in a separate working
folder outside this code repository and are not tracked here.

---

## Execution Order (v2, current)

```r
# Step 1. Model qualification (parent/metabolite Css separated)
source("R/qualification/qualification_v2.R")

# Step 2. Monte Carlo simulation (reuses the frozen n=10,000/drug cohort;
#          test with N <- 100 first, then N <- 10000 for the real run —
#          see data/pk_parameters/pk_parameters_v2.R and
#          R/simulation/mc_run_v2.R header comments)
source("R/simulation/mc_run_v2.R")

# Step 3. Range-sensitivity analyses (fm,CYP2D6 sweep + metabolite potency sweep)
source("R/sensitivity/sensitivity_v2.R")

# Step 4. Summary tables (Table 3-6, S1, B1, Dunn tests)
source("R/analysis/make_tables_v2.R")

# Step 5. Figures
source("R/figures/figure1_discontinuation_v2.R")
source("R/figures/figure2_venlafaxine_IR_XR_v2.R")
source("R/figures/figure3_mc_t50_violin_v2.R")
source("R/figures/figure4_mc_dCdt_ridge_v2.R")
source("R/figures/figure5_mc_scatter_v2.R")
```

Full step-by-step instructions, expected runtimes, and troubleshooting for
each v2 script are in `04_resubmission/RUN_ORDER.md`.

---

## Drug & PK Parameters (v2)

| Drug | Class | Model | Parent t½ (label) | F |
|---|---|---|---|---|
| Sertraline | SSRI, IR | 1-cmt | 26 h | 1.0 |
| Fluoxetine + norfluoxetine | SSRI, IR | 1-cmt + metabolite | 153 h / 223 h | 1.0 |
| Paroxetine | SSRI, IR | 1-cmt | 21 h | 1.0 |
| Venlafaxine IR + ODV | SNRI, IR | 1-cmt + metabolite | 5 h / 11 h | 1.0 |
| Venlafaxine XR + ODV | SNRI, XR | 1-cmt + metabolite | 5 h / 11 h | 1.0 |

Bioavailability is unified to F = 1.0 for all five regimens (D3); see
`data/pk_parameters/pk_parameters_v2.R` and manuscript Methods 2.2/2.3 for
full parameter values, units, and literature/author-estimate provenance.
Note: the specific primary source for the fluoxetine CL/Vd per-kg figures
and the paroxetine Vd per-kg figure could not be reconfirmed during the
v2 revision (flagged in the manuscript, Methods 2.2 and Limitations).

---

## Key PK Metrics (v2 definitions)

| Metric | Definition |
|---|---|
| C0 | Instantaneous plasma concentration at the moment of discontinuation (D1 — not an interval-averaged Css) |
| t50% C0 | Time to reach 50% of C0 (linear interpolation); reached/censored/calc_failed status coded per patient |
| t25% C0 | Time to reach 25% of C0 |
| \|Rel. dC/dt\| @24h | \|(C0 − C24h)\| / (24 × C0) |
| AUC deficit (absolute) | C0 × (observation window) − area under the post-discontinuation curve |
| AUC deficit (normalised) | 100 × absolute AUC deficit / (C0 × observation hours), % — directly comparable across drugs regardless of absolute concentration scale |

---

## Dependencies

```r
install.packages(c(
  "rxode2",     # ODE solver
  "tidyverse",  # Data manipulation
  "patchwork",  # Figure assembly
  "ggridges",   # Ridge plots
  "ggrepel",    # Label repulsion
  "rstatix"     # Statistical tests
))
```

---

## Reproducibility

- Random seed: `set.seed(2025)` applied in all stochastic steps
- Monte Carlo: n = 10,000 simulated patients per drug (50,000 total across
  5 regimens), frozen cohort reused for both v1 and v2 analyses (sampling
  itself did not change — only the C0/F/window definitions applied to it)
- Stability confirmed: CV < 5% at n=1,000; CV < 2% at n=5,000
- R version: 4.3.3 | rxode2: 5.0.2 | ggplot2: 4.0.3

---

## Notes

- v2 uses the instantaneous C0 at discontinuation (D1), not interval-averaged Css, as the reference concentration for all time-to-threshold metrics
- Bioavailability F = 1.0 unified across all regimens (D3); legacy per-drug F values from v1 are documented but no longer applied
- Active metabolite contribution uses an unweighted sum as the base case; a potency-weighted sensitivity analysis is reported separately (manuscript Section 3.9)
- CYP2D6 variability applied to fluoxetine, paroxetine, venlafaxine IR/XR via a uniform drug-independent clearance multiplier (explicitly framed as a simplified sensitivity analysis, not a drug-specific pharmacogenomic model)
- No patients are excluded from the Monte Carlo analysis; all are status-coded as reached/censored/calc_failed (replaces the v1 `t50% < 1h` exclusion rule)
- All v1 scripts and outputs are retained unmodified for audit purposes; `_v2` files are authoritative for the resubmitted manuscript

---

## AI-assisted development disclosure

Portions of this codebase (simulation scripts, statistical analysis, and
documentation) were developed with the assistance of Claude (Anthropic).
All AI-assisted code was reviewed, executed, and verified by the authors
against its console output and the manuscript's reported numbers. See the
manuscript's Declarations section for the full AI-use disclosure.

---

*Last updated: 2026-10-26 (v2 resubmission)*
