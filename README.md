# Pharmacokinetic Simulation of Post-Discontinuation Plasma Concentration Decline Across SSRI and SNRI Regimens: An Exploratory Monte Carlo Framework

**Authors:** Nari Kim, Dukhyun Kim
**Institution:** Kingston University London, BSc Pharmacology
**Target Journal:** Translational and Clinical Pharmacology (TCP), manuscript TCP-26-013
**Status:** Revision 1 resubmission (v3 re-analysis) — tag `tcp-26-013-revision1`

---

## Research Overview

This repository contains the computational pharmacokinetic (PK) simulation
framework used to compare post-discontinuation plasma concentration decline
across five SSRI/SNRI regimens (sertraline, fluoxetine + norfluoxetine,
paroxetine, venlafaxine IR, venlafaxine XR) under abrupt discontinuation,
including a paroxetine CYP2D6 poor-metabolizer clearance ratio as a model
input.

### Versions in this repository

- **v1** — code and results as they stood at the original submission
  (2026-09).
- **v2** — an internal revision that was never submitted; kept for audit
  only (see `outputs/tables_v2/s1_v2_old_vs_new_paired.csv` for a
  patient-level comparison against v1).
- **v3** — the revision-1 resubmission. **This is the version the
  manuscript cites.** All `_v3` files (and `_v3_parV835`, a full-cohort
  rerun used only for the paroxetine V/F coefficient-of-variation
  sensitivity analysis, manuscript Section 3.9) are current; `_v2` and
  unsuffixed (`v1`) files are retained for audit only and should not be
  cited as current results.

### Key Findings (v3, as reported in the manuscript Abstract)

- Median time to 50% of the discontinuation concentration (t₅₀% C₀) was
  16.0 h for venlafaxine XR, 16.5 h for IR, 21.9 h for paroxetine, 27.0 h
  for sertraline and 323 h for fluoxetine.
- For the three regimens with an active metabolite, t₅₀% C₀ was 2.09 to
  3.04 times the implemented parent half-life.
- The 7-day normalised AUC deficit was 15.6% for fluoxetine and 77.1% to
  86.8% for the other regimens.
- Paroxetine poor metabolizers had a longer median t₅₀% C₀ than normal
  metabolizers (34.4 vs 21.0 h; rank-biserial r = +0.458).
- In a published placebo-substitution study, the fraction remaining after
  5 days of placebo was close to the simulated value for fluoxetine and
  higher than simulated for paroxetine.

The model describes plasma concentrations only and was not tested against
clinical discontinuation symptom outcomes. See the manuscript Abstract and
Results for full detail; this README does not restate numbers beyond what
the manuscript itself reports.

---

## Repository Structure

```
.
├── R/
│   ├── models/                     # v1 deterministic PK models (unchanged)
│   ├── qualification/
│   │   ├── model_qualification.R   # v1 (superseded)
│   │   ├── qualification_v2.R      # v2 (not submitted; audit only)
│   │   └── qualification_v3.R      # v3 — current (Table 2a/2b)
│   ├── simulation/
│   │   ├── mc_params.R / mc_run_single.R     # v1 (superseded)
│   │   ├── mc_run_v2.R                        # v2 (audit only)
│   │   ├── mc_params_v3.R / mc_run_v3.R       # v3 — current (n=10,000/drug)
│   │   └── mc_params_v3_parV835.R / mc_run_v3_parV835.R
│   │       # v3, full cohort rerun with paroxetine V/F CV = 83.5% instead
│   │       # of the primary 35% — used only for the Section 3.9 V/F CV
│   │       # sensitivity analysis, same virtual patients as the primary run
│   ├── sensitivity/
│   │   ├── sensitivity_analysis.R / sensitivity_v2.R   # v1/v2 (superseded/audit)
│   │   └── sensitivity_v3.R         # v3 — current (poor-metabolizer ratio
│   │                                 #  and metabolite potency weight sweeps)
│   ├── analysis/
│   │   ├── stability_check.R        # v1 (superseded)
│   │   ├── stability_v3.R           # v3 — current (sample-size stability,
│   │   │                            #  Supplementary Table S5)
│   │   ├── phase5_statistics.R / make_tables_v2.R   # v1/v2 (superseded/audit)
│   │   └── tables_v3.R              # v3 — current (Tables 2a/2b/3/4,
│   │                                 #  Supplementary Tables S1-S4)
│   └── figures/
│       ├── figure1-5_*.R / *_v2.R   # v1/v2 (superseded/audit)
│       └── figure1-5_v3.R           # v3 — current, generates Figures 1-5
│                                     #  reported in the manuscript
│
├── data/
│   └── pk_parameters/
│       ├── pk_parameters.R / pk_parameters_v2.R   # v1/v2 (superseded/audit)
│       ├── pk_parameters_v3.R         # v3 — current (Table 1)
│       └── pk_parameters_v3_parV835.R # v3, V/F CV = 83.5% variant (Section 3.9)
│
├── outputs/
│   ├── figures/        # v1/v2 PNGs kept for audit; v3 figures ship alongside
│   │                    # this release as outputs/figures/*_v3.png
│   ├── qualification/ , qualification_v2/ , qualification_v3/
│   ├── mc/ , mc_v2/ , mc_v3/ (incl. *_parV835 variant)
│   ├── mc_sensitivity/ , mc_sensitivity_v2/ , mc_sensitivity_v3/
│   ├── stability_v3/
│   └── tables/ , tables_v2/ , tables_v3/
```

Manuscript drafts, the response letter, cover letter, and the full working
record (handoff notes, independent-calculation cross-checks, literature PDFs)
live in a separate working folder outside this code repository and are not
tracked here.

---

## Execution Order (v3, current)

```r
# Step 1. PK parameters (Table 1)
source("data/pk_parameters/pk_parameters_v3.R")

# Step 2. Model qualification (Table 2a/2b)
source("R/qualification/qualification_v3.R")

# Step 3. Monte Carlo simulation (n = 10,000 virtual patients per drug,
#          50,000 total; primary cohort)
source("R/simulation/mc_params_v3.R")
source("R/simulation/mc_run_v3.R")

# Step 3b. V/F CV = 83.5% full-cohort rerun (same virtual patients),
#           used only for the Section 3.9 V/F CV sensitivity analysis
source("data/pk_parameters/pk_parameters_v3_parV835.R")
source("R/simulation/mc_params_v3_parV835.R")
source("R/simulation/mc_run_v3_parV835.R")

# Step 4. Range-sensitivity analyses (poor-metabolizer clearance ratio,
#          metabolite potency weight; n = 2,000 per condition)
source("R/sensitivity/sensitivity_v3.R")

# Step 5. Sample-size stability analysis (Supplementary Table S5)
source("R/analysis/stability_v3.R")

# Step 6. Summary tables (Table 2a/2b/3/4; Supplementary Tables S1-S4)
source("R/analysis/tables_v3.R")

# Step 7. Figures
source("R/figures/figure1-5_v3.R")
```

A fixed seed is set in each script; see script headers for the exact value
used in each step.

**Note for reproduction on another machine:** every `R/*_v3*.R` and
`data/pk_parameters/*_v3*.R` script hardcodes
`BASE <- "~/nari-research/pkpd-antidepressant-sim"` as the repository root.
Change this line to your own local path before running.

---

## Dependencies

```r
install.packages(c(
  "rxode2",     # ODE solver (version 5.1.8 used for v3; see manuscript ref. [18])
  "tidyverse",  # Data manipulation
  "patchwork",  # Figure assembly
  "ggridges",   # Ridge plots
  "ggrepel",    # Label repulsion
  "rstatix"     # Statistical tests
))
```

R and package versions used for the v3 run were not separately recorded in
a saved `sessionInfo()` output; the v2 run (same machine) used R 4.3.3.
rxode2 version 5.1.8 is cited in the manuscript (reference [18]).

---

## Parameter Sourcing

Full parameter values, units and literature/author-assumption provenance
for v3 are in the manuscript, Methods Section 2.2 and Table 1 — this
README does not restate them separately so that the two cannot drift out
of sync. For a machine-readable copy of the same parameter set, see
`data/pk_parameters/pk_parameters_v3.R`.

---

## Reproducibility

- Monte Carlo: n = 10,000 simulated patients per drug (50,000 total across
  5 regimens) for the v3 primary analysis and for the V/F CV = 83.5%
  sensitivity rerun (same virtual-patient cohort, recovered by inverting
  each patient's earlier standard-normal draw — manuscript Methods
  Section 2.4.2); n = 2,000 per condition for the two other
  range-sensitivity analyses (manuscript Section 3.9).
- Sample-size stability: CV of the median t₅₀% C₀ across resampling seeds
  fell to ≤1.5% at n = 10,000 for every drug (Supplementary Table S5).
- See script headers for the exact seed used in each step.

---

## AI-assisted development disclosure

Portions of this codebase (simulation scripts, statistical analysis, and
documentation) were developed with the assistance of Claude (Anthropic).
All AI-assisted code was reviewed, executed, and verified by the authors
against its console output and the manuscript's reported numbers. See the
manuscript's Acknowledgments section for the full AI-use disclosure.

---

*Last updated: 2026-10-30 (v3 revision-1 resubmission)*
