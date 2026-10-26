# ============================================================
# PK Parameters v2 — Antidepressant Discontinuation Risk Study
# Author   : Nari Kim (코드 작성 보조: Claude)
# Based on : data/pk_parameters/pk_parameters.R (v4, 2025-05)
#
# v2에서 바뀐 것 (원본 파일은 그대로 두고 새 파일로 만듦 — 덕현 인계
# HANDOFF_TO_NARI.md 6-3, RESUBMISSION_PLAN.md 2단계 D2·D3 결정 반영)
#
#   D3 (확정): 모든 약물의 F를 1.0으로 통일.
#     원래 F<1 값(fluoxetine 0.70, paroxetine 0.50, venlafaxine 0.45)은
#     CL/Vd가 이미 라벨의 "apparent clearance(CL/F)"·"apparent volume(V/F)"
#     — 즉 F로 나뉜 값 — 라는 추정 하에 용량에 F를 한 번 더 곱해서
#     생체이용률을 이중 적용하고 있었을 가능성이 높음(CHECK_RESULTS_1026.md
#     6절). F=1로 두면 투여량이 전체 공칭 용량(dose_mg)으로 들어가고,
#     CL·Vd는 "이미 F로 나뉜 apparent 값"으로 그대로 해석됨 — 기존 숫자를
#     바꾸는 게 아니라 중복 적용만 없앤 것.
#
#   D2 (확정): CL·Vd·Ka 등 숫자 값 자체는 전부 유지. 다만:
#     - 각 파라미터에 basis 태그 추가: "label"(라벨 직접 보고값) /
#       "literature"(타 문헌) / "derived"(반감기 등에서 역산) /
#       "unsourced_pending"(출처 미확인, SOURCE_CHECKS.md PENDING과 동일)
#     - sertraline Vd의 "20 L/kg, Monfort 2024 BJCP" 표기를 제거함.
#       덕현 쪽 점검에서 Monfort 2024 원문 Table 2는 V/F=2,320 L,
#       CL/F=69.6–119 L/h로 보고돼 있고 "20 L/kg"이라는 값 자체가 없다는
#       점이 확인됨(HANDOFF_TO_NARI.md 3-6). 숫자(1400 L)는 그대로 두고
#       출처만 "unsourced_pending"으로 내림 — 실제 출처를 못 찾으면 원고에
#       저자 가정으로 명시해야 함.
#     - CYP2D6 배수(0.30/1.00/1.60)에도 출처가 없다는 점을 명시적으로
#       주석에 남김(HANDOFF_TO_NARI.md 3-6) — 숫자는 유지.
#     - venlafaxine 구현 반감기(Vd/CL 역산 ≈ 4.0 h)가 Table 1에 적힌
#       라벨 값 5 h와 다르다는 점을 주석으로 남김 — Table 1 각주 작성용.
#
# 이 파일은 mc_run_v2.R, qualification_v2.R, sensitivity_v2.R이
# 공통으로 source() 해서 쓴다. F 값을 스크립트마다 따로 적지 않기 위함
# (HANDOFF_TO_NARI.md 8-1: "F 값은 pk_parameters.R 한 곳에서만 바꾸고,
# 모든 스크립트가 그 파일을 읽는지 확인" — sensitivity_analysis.R이
# F를 직접 하드코딩했던 문제의 재발 방지).
# ============================================================

pk_params <- list(

  # ── Sertraline ─────────────────────────────────────────────
  # CL = 0.693 × 1400 / 26 = 37.3 L/h (FDA label t½=26h 역산)
  sertraline = list(
    drug        = "Sertraline",
    class       = "SSRI",
    formulation = "IR",
    Ka          = 0.50,    # 1/h | Tmax 4.5–8.4h 역산 | FDA label
    Ka_basis    = "derived",
    Vd          = 1400,    # L   | 숫자는 v4와 동일하게 유지
    Vd_basis    = "unsourced_pending",  # "20 L/kg, Monfort 2024"는 원문과 불일치 확인됨 — 출처 재확인 필요
    CL          = 37.3,    # L/h | 0.693 × 1400 / 26 | FDA label t½ 역산
    CL_basis    = "derived",
    F           = 1.0,     # v2: D3 — 전 약물 F=1 통일 (원래도 1.0, 변경 없음)
    F_basis     = "apparent (F 미확립, CL/F·V/F로 해석)",
    t_half      = 26,      # h   | FDA label
    CV_CL       = 0.40,
    CV_Vd       = 0.35,
    model       = "1-compartment",
    source      = "FDA label (Zoloft, t_half); Vd 출처는 재확인 필요(PENDING, Monfort 2024 원문과 불일치)"
  ),

  # ── Fluoxetine ─────────────────────────────────────────────
  fluoxetine = list(
    drug        = "Fluoxetine",
    class       = "SSRI",
    formulation = "IR",
    Ka          = 0.72,    # 1/h | Tmax 7.1h 역산 | FDA label 6-8h 기준
    Ka_basis    = "derived",
    Vd          = 2310,    # L   | 33 L/kg × 70kg | NZ Medicines
    Vd_basis    = "literature_not_in_reflist",  # NZ Medicines — 참고문헌 목록에 없음, 추가 또는 대체 필요
    CL          = 10.5,    # L/h | 0.15 L/h/kg × 70kg | NZ Medicines
    CL_basis    = "literature_not_in_reflist",
    F           = 1.0,     # v2: D3 — 원래 0.70 → 1.0 (CL/F·V/F 이중적용 제거)
    F_basis     = "apparent (CL/F·V/F로 해석, v4의 F=0.70 중복 적용 제거)",
    t_half      = 153,     # h   | 6.4일 (CL=10.5 역산값, 문헌범위 4-16일 내)
    CV_CL       = 0.874,
    CV_Vd       = 0.50,
    model       = "1-compartment + metabolite (단순화, terminal elimination 집중)",
    norfluoxetine = list(
      Ka_met    = 0.30,
      Vd_met    = 2310,
      Vd_met_basis = "literature_not_in_reflist",
      CL_met    = 7.18,    # L/h | 0.693 × 2310 / 223 | FDA label t½=223h 역산
      CL_met_basis = "derived",
      t_half    = 223,     # h   | FDA label 반복투여 9.3일
      CV_CL_met = 0.564,
      Fm        = 0.72,
      Fm_basis  = "unsourced_pending"  # 저자 가정, SOURCE_CHECKS.md PENDING
    ),
    source      = "NZ Medicines(CL/Vd, 참고문헌 목록에 없음); FDA label(Ka/t½_nfx); Pauchard 2011(Tmax기준)"
  ),

  # ── Paroxetine IR ──────────────────────────────────────────
  paroxetine = list(
    drug        = "Paroxetine",
    class       = "SSRI",
    formulation = "IR",
    Ka          = 0.580,   # 1/h | Tmax 5.2h 역산 | Li 2022 Tmax 기준
    Ka_basis    = "derived",
    Vd          = 609,     # L   | 8.7 L/kg × 70kg | DrugBank/medicine.com
    Vd_basis    = "literature_not_in_reflist",  # DrugBank — 참고문헌 목록에 없음
    CL          = 20.1,    # L/h | 0.693 × 609 / 21 | FDA label t½ 역산
    CL_basis    = "derived",
    F           = 1.0,     # v2: D3 — 원래 0.50 → 1.0 (CL/F·V/F 이중적용 제거)
    F_basis     = "apparent (CL/F·V/F로 해석, v4의 F=0.50 중복 적용 제거, Kaye 1989도 참고문헌 목록에 없음)",
    t_half      = 21,      # h   | FDA label (범위 7–65h)
    CV_CL       = 0.60,    # IIV 높음 (비선형 CYP2D6)
    CV_Vd       = 0.35,
    model       = "1-compartment",
    CL_note     = "CL = FDA label t½ 역산. Li 2022 CL/F=21.2 미적용 (V/F=8850 왜곡)",
    source      = "FDA label (Paxil, t_half); DrugBank(Vd, 참고문헌 목록에 없음); Li 2022(Ka/Tmax만)"
  ),

  # ── Venlafaxine IR ─────────────────────────────────────────
  # 참고: Vd/CL로 역산한 구현 반감기 = 0.693×525/91 ≈ 4.0h.
  #       Table 1에 적힌 라벨 값 5h와 다름 — Vd, CL을 각각 독립적으로
  #       반올림/추정했기 때문에 생기는 내부 불일치. 숫자는 유지하고
  #       Table 1 각주에 이 차이를 명시해야 함(HANDOFF_TO_NARI.md 3-6).
  venlafaxine_IR = list(
    drug        = "Venlafaxine",
    class       = "SNRI",
    formulation = "IR",
    Ka          = 1.50,    # 1/h | Tmax ~2h | FDA label
    Ka_basis    = "label",
    Vd          = 525,     # L   | 7.5 L/kg × 70kg | FDA label
    Vd_basis    = "label",
    CL          = 91,      # L/h | 1.3 L/h/kg × 70kg | FDA label
    CL_basis    = "label",
    F           = 1.0,     # v2: D3 — 원래 0.45 → 1.0 (CL/F·V/F 이중적용 제거)
    F_basis     = "apparent (CL/F·V/F로 해석, v4의 F=0.45 중복 적용 제거). 라벨: 흡수율 92%, 절대생체이용률 45% → CL 1.3 L/h/kg는 CL/F로 해석됨(HANDOFF_TO_NARI.md 3-4, 추론)",
    t_half      = 5,       # h   | FDA label (구현 반감기 Vd/CL≈4.0h와 다름, 각주 필요)
    CV_CL       = 0.46,
    CV_Vd       = 0.49,
    model       = "1-compartment + metabolite (ODV)",
    ODV = list(
      Vd_odv    = 399,
      CL_odv    = 28,
      t_half    = 11,
      CV_CL_odv = 0.50,
      Fm        = 0.80,
      Fm_basis  = "unsourced_pending"
    ),
    source      = "FDA label (Effexor XR 2012)"
  ),

  # ── Venlafaxine XR ─────────────────────────────────────────
  venlafaxine_XR = list(
    drug        = "Venlafaxine",
    class       = "SNRI",
    formulation = "XR",
    Ka          = 0.25,    # 1/h | Tmax ~5.5h | FDA label
    Ka_basis    = "label",
    Vd          = 525,
    Vd_basis    = "label",
    CL          = 91,
    CL_basis    = "label",
    F           = 1.0,     # v2: D3 — 원래 0.45 → 1.0
    F_basis     = "apparent (CL/F·V/F로 해석, v4의 F=0.45 중복 적용 제거)",
    t_half      = 5,
    CV_CL       = 0.46,
    CV_Vd       = 0.49,
    model       = "1-compartment + metabolite (ODV)",
    ODV = list(
      Vd_odv    = 399,
      CL_odv    = 28,
      t_half    = 11,
      CV_CL_odv = 0.50,
      Fm        = 0.80,
      Fm_basis  = "unsourced_pending"
    ),
    source      = "FDA label (Effexor XR 2012)"
  )
)

# ── CYP2D6 genotype 보정 ───────────────────────────────────
# 배수(0.30/1.00/1.60) 자체의 출처는 확인되지 않음(HANDOFF_TO_NARI.md
# 3-6) — 숫자는 v4와 동일하게 유지하고 Limitation에 명시 필요.
# "rapid" 22%가 Gaedigk 2017의 phenotypic/genotype-predicted 중 어느
# 쪽인지도 확인 필요(HANDOFF_TO_NARI.md 8-3, B6).
cyp2d6 <- list(
  poor_metabolizer   = list(freq = 0.08, CL_factor = 0.30, freq_basis = "Gaedigk 2017 — 정확한 범주(phenotypic vs genotype-predicted) 확인 필요", CL_factor_basis = "unsourced_pending"),
  normal_metabolizer = list(freq = 0.70, CL_factor = 1.00, freq_basis = "Gaedigk 2017 — 정확한 범주 확인 필요", CL_factor_basis = "reference (배수 1.0)"),
  rapid_metabolizer  = list(freq = 0.22, CL_factor = 1.60, freq_basis = "Gaedigk 2017 — 정확한 범주 확인 필요", CL_factor_basis = "unsourced_pending")
)

# ── 검증 요약 ──────────────────────────────────────────────
cat("✅ PK 파라미터 v2 로드 완료 (v4 기반, D2 출처 표기 정정 + D3 F=1 통일)\n")
cat("약물:", paste(names(pk_params), collapse=" | "), "\n\n")
cat("[v2 변경사항 — v4 대비]\n")
cat("  F: fluoxetine 0.70→1.0, paroxetine 0.50→1.0, venlafaxine IR/XR 0.45→1.0 (D3, 이중적용 제거)\n")
cat("  sertraline Vd 출처 표기: 'Monfort 2024 BJCP'(원문과 불일치) → 'unsourced_pending' (D2, 숫자 1400L은 유지)\n")
cat("  각 파라미터에 basis 태그 추가 (label/literature_not_in_reflist/derived/unsourced_pending)\n")
cat("  CL, Vd, CV, Ka 등 숫자 값 자체는 v4와 전부 동일 (D2: 값 유지)\n")
