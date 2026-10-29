# ============================================================
# PK Parameters v3 -- parV835 변형 (게이트 전용, 본 분석 아님)
# Author   : Nari Kim (code assist: Claude)
# Based on : data/pk_parameters/pk_parameters_v3.R
#
# 이 파일은 pk_parameters_v3.R을 그대로 복사한 뒤 paroxetine CV_Vd 한
# 값만 0.35 -> 0.835로 바꾼 변형이다 (HANDOFF_v3_ADDENDUM_1028 3-4 7번,
# self_check_v3.py의 --expected-params expected_params_v3_parV835.csv /
# --expected-mc expected_mc_v3_parV835.csv against용). 본 분석(pk_parameters_v3.R,
# CV_Vd=0.35)은 건드리지 않음 -- 이 파일과 mc_params_v3_parV835.R,
# mc_run_v3_parV835.R은 게이트 확인이 끝나면 더 쓰지 않아도 되는 1회용
# 보조 스크립트다 (HANDOFF_v3_ADDENDUM_1028의 "구조를 임의로 바꾸지 말고
# 알려 주세요" 원칙에 따라, 이미 PASS한 메인 파이프라인 스크립트/산출물은
# 그대로 두고 별도 파일로만 비교용 변형을 만듦).
#
# v3에서 바뀐 것 (HANDOFF_TO_NARI_v3.md 5절 Table 1, D2'.D5.D6.D7.D9 결정 반영)
#
#   D2' (확정, v2의 D2를 대체): 추정값을 라벨 재보정값으로 교체.
#     v2는 "숫자는 유지, 출처 표기만 정정"이었으나, v3는 sertraline,
#     paroxetine, venlafaxine의 CL/F, V/F 숫자 자체를 라벨 1차 자료에
#     맞춰 다시 계산함. fluoxetine 값은 유지하고 출처만 교체.
#
#   D5 (확정): CYP2D6 보정은 paroxetine에만 적용. PM x 0.58
#     (Sindrup 1992, Clin Pharmacol Ther 51:278-87, Table II, 정상상태
#     30 mg/day 총경구 청소율 중앙값 비 EM 36 / PM 21 L/h = 0.58).
#     RM 범주는 삭제(근거 수준 부족 - evidence/CYP2D6_EVIDENCE_1027.md
#     3절: "RM(1.60, 빈도 22%) ... 빼는 것을 권합니다"). fluoxetine,
#     venlafaxine은 CYP2D6를 모델링하지 않음 - 라벨/CPIC 2023 보충자료
#     근거는 "대사체 비는 바뀌지만 합계(parent+metabolite)는 크게
#     바뀌지 않는다"이고, 단순 CL 배수 구조로는 이 구조를 재현할 수
#     없어서 모델링하지 않음("phenotype과 무관"이라고 쓰면 안 됨).
#
#   D6 (확정): paroxetine 임상 투여 설계를 30 mg/day로 올림
#     (censoring 감소 목적, HANDOFF_v3 6절).
#
#   D7 (확정): venlafaxine parent, ODV CL/F를 Effexor XR 150 mg QD
#     라벨 정상상태 Cmax(parent 150, ODV 260 ng/mL)에 맞춰
#     calibration/ven_calibration_v3.py(brentq)로 보정. V/F는 라벨값
#     고정(7.5 L/kg parent, 5.7 L/kg ODV). venlafaxine IR 75mg BID의
#     라벨 Cmax(225/290)는 보정에 쓰지 않고 독립 점검용으로만 보고.
#
#   D9 (확정): CV 출처 재조정 - Monfort 2024(sertraline CL),
#     Li 2022(paroxetine), 라벨 SD/평균(venlafaxine, ODV) 재확인.
#     Sagahon-Azua 2021 Table 2의 fluoxetine/norfluoxetine CV는
#     그 논문의 Bayesian 개별추정에 쓰인 "사전값(prior)"이지 최종
#     집단 추정치가 아니라는 점을 주석과 원고 모두에 명시해야 함
#     (HANDOFF_v3 5절 "주의할 출처" 2번째 항목).
#
#   Gaedigk 2017 PM 빈도 재확인(HANDOFF_v3 5절 지시 이행, 2026-10-28
#   WebFetch로 원문 확인):
#     Gaedigk, Sangkuhl, Whirl-Carrillo et al. "Prediction of CYP2D6
#     phenotype from genotype across world populations." Genet Med.
#     2017;19:69-76. Table 1 / Supplementary Table S2(유전형 예측치)의
#     유럽인 평균 PM 빈도는 5.4%이고, 개별 하위군 중 Ashkenazi Jewish만
#     6%로 8%에 가장 가까움 - 8%와 정확히 일치하는 유전형 예측 하위군은
#     없음. 반면 Discussion 문단은 유럽인의 "관찰된 표현형(phenotypic)
#     PM 빈도(pPM)"를 8.45%로 보고하며, 이 값이 기존 8%와 사실상 일치함.
#     따라서 8%는 유지하되, 출처를 "Table 1(유전형 예측, 5.4%)"이 아니라
#     "Discussion(표현형 관찰, pPM 8.45%, Europeans)"으로 정확히 교정함.
#     원고 각주에도 이 구분(표현형 vs 유전형 예측)을 명시할 것.
#
#   삭제한 것(HANDOFF_v3 5절 "삭제할 것"):
#     - norfluoxetine Ka_met 0.30 (모델 식에 대사체 Ka가 없음, B4)
#     - sertraline Vd 1,400 L (라벨 재보정값 3,054 L로 교체)
#     - CYP2D6 rapid_metabolizer(RM) 범주 전체
#     - DrugBank 인용, Medsafe 인용(fluoxetine 출처를 Christchurch
#       Medicines Information Service(1999)로 교체)
#
# 이 파일은 mc_run_v3.R, qualification_v3.R, sensitivity_v3.R,
# stability_v3.R이 공통으로 source() 해서 쓴다.
# ============================================================

pk_params <- list(

  # -- Sertraline --------------------------------------------
  # CL/F = 200 mg / 2.570 mg.h/L = 77.8 L/h
  #   (Zoloft 2009 FDA label, 200 mg/day 30일 정상상태 평균 AUC0-24
  #   2,570 ng.h/mL)
  # V/F = CL/F x t1/2 / ln2 = 77.8 x 27.2 / 0.693 ~= 3,054 L
  #   (같은 성인군 평균 t1/2 27.2h에서 유도 - Zoloft 2009 label)
  sertraline = list(
    drug        = "Sertraline",
    class       = "SSRI",
    formulation = "IR",
    Ka          = 0.50,      # 1/h | 라벨 Tmax에서 유도 (v2와 값 동일, 유지)
    Ka_basis    = "derived",
    Vd          = 3054,      # L (V/F) | CL/F x t_half / ln2 | Zoloft 2009 label
    Vd_basis    = "label_derived",
    CL          = 77.8,      # L/h (CL/F) | 200/2.570 | Zoloft 2009 label AUC0-24
    CL_basis    = "label_derived",
    F           = 1.0,
    F_basis     = "apparent (CL/F, V/F로 해석)",
    t_half      = 27.2,      # h | V/F, CL/F에서 유도 (위 두 값과 내적 일치)
    CV_CL       = 0.415,     # Monfort 2024 [7] Table 2
    CV_Vd       = 0.35,      # 저자 가정 (Monfort는 V/F IIV를 추정하지 않음)
    CV_Vd_basis = "unsourced_pending (author assumption)",
    model       = "1-compartment",
    source      = "Zoloft 2009 FDA label (CL/F, V/F, Ka); Monfort 2024 [7] Table 2 (CV_CL)"
  ),

  # -- Fluoxetine ----------------------------------------------
  fluoxetine = list(
    drug        = "Fluoxetine",
    class       = "SSRI",
    formulation = "IR",
    Ka          = 0.72,      # 1/h | v2와 값 동일, 유지
    Ka_basis    = "derived",
    Vd          = 2310,      # L (V/F) | 33 L/kg x 70kg | 값 유지, 출처만 교체
    Vd_basis    = "literature_secondary",
    CL          = 10.5,      # L/h (CL/F) | 0.15 L/h/kg x 70kg | 값 유지, 출처만 교체
    CL_basis    = "literature_secondary",
    F           = 1.0,
    F_basis     = "apparent (CL/F, V/F로 해석)",
    t_half      = 152.5,     # h | 값 유지
    CV_CL       = 0.874,     # Sagahon-Azua 2021 [6] Table 2 -- 주의: 그 논문의
                              # Bayesian 개별추정 "사전값"이지 최종 집단 추정치가
                              # 아님 (HANDOFF_v3 5절 주의사항 2번째). 2026-10-27
                              # 재검증으로 원출처가 Panchaud 2011 (Clin Pharmacol
                              # Ther 2011;89(6):830-6, DOI 10.1038/clpt.2011.23)임을
                              # 확인함 -- Sagahon-Azua 2021이 Bayesian 사전값으로
                              # 차용한 것. evidence/CLAUDE_3PAPERS_REVERIFICATION_1027.md
    CV_CL_note  = "prior value originally from Panchaud et al. 2011 (Clin Pharmacol Ther 89(6):830-6), adopted as a Bayesian prior by Sagahon-Azua et al. 2021 Table 2 -- NOT a population estimate from either paper's own fitted cohort; cite Panchaud 2011 accordingly in manuscript",
    CV_Vd       = 0.50,      # 저자 가정
    CV_Vd_basis = "unsourced_pending (author assumption)",
    model       = "1-compartment + metabolite (단순화, terminal elimination 집중)",
    norfluoxetine = list(
      # Ka_met 삭제됨 (모델 식에 대사체 Ka가 없음, HANDOFF_v3 5절 "삭제할 것" B4)
      Vd_met       = 2310,
      Vd_met_basis = "값 유지 (literature_secondary)",
      CL_met       = 7.18,   # L/h | 0.693 x 2310 / 223 | 값 유지
      CL_met_basis = "derived",
      t_half       = 223,    # h | 값 유지
      CV_CL_met    = 0.564,  # Sagahon-Azua 2021 Table 2 -- 위와 동일한 주의
                              # 사항 (prior, not final estimate; 원출처 Panchaud
                              # 2011 -- 위 CV_CL_note 및 evidence 파일 참고)
      Fm           = 0.72,
      Fm_basis     = "unsourced_pending (author assumption)"
    ),
    source        = "Christchurch Medicines Information Service \"FLUOXETINE\" (1999) -- V 33 L/kg, CL 0.15 L/h/kg (2차 자료, compiled TDM summary; 원출처로 보이는 Benfield 1986 Drugs 32:481-508 원문 미확보); Sagahon-Azua 2021 [6] Table 2 (CV, prior value, originally from Panchaud 2011 Clin Pharmacol Ther 89(6):830-6)",
    source_caveat = "Christchurch 자료에는 F=0.9가 별도로 적혀 있어 33 L/kg, 0.15 L/h/kg이 전신(systemic) 값일 가능성이 있음 -- 원고에 'we interpreted these as apparent values; if they were systemic values with F=0.9, CL/F would be about 11% higher' 취지를 명시 (HANDOFF_v3 5절 주의사항 1번째)"
  ),

  # -- Paroxetine IR -------------------------------------------
  paroxetine = list(
    drug        = "Paroxetine",
    class       = "SSRI",
    formulation = "IR",
    Ka          = 0.39,      # 1/h | Tmax 5.2h에서 유도 | Paxil 2017 label
    Ka_basis    = "label_derived",
    Vd          = 813,       # L (V/F) | Cmax기준값, Cmin기준값의 기하평균 | Paxil 2017 label
    Vd_basis    = "label_derived",
    CL          = 26.8,      # L/h (CL/F) | 위와 동일 방법
    CL_basis    = "label_derived",
    F           = 1.0,
    F_basis     = "apparent (CL/F, V/F로 해석)",
    t_half      = 21.0,      # h | Paxil 2017 label 30mg/day 30일 정상상태
    CV_CL       = 0.519,     # Li 2022 [8] Table 2 -- 대상: 중국 정신질환
                              # 환자 184명 (15-90세); CYP2D6 공변량 미포함 시
                              # IIV에 phenotype 변이가 이미 포함되어 PM
                              # 배수와 겹칠 수 있음 -- Limitations(B10)에
                              # 명시 필요
    CV_Vd       = 0.835,     # parV835 변형 -- 본 분석 값(0.35)이 아니라 Li 2022
                              # Table 2의 V/F IIV 원값 83.5%를 그대로 씀.
                              # self_check_v3.py 게이트(expected_params_v3_parV835.csv
                              # / expected_mc_v3_parV835.csv) 대조 전용. 본 분석
                              # 결론(CV_Vd=0.35, pk_parameters_v3.R)에는 영향 없음.
    CV_Vd_basis = "Li 2022 [8] Table 2 원값, 83.5% (parV835 게이트 대조 전용 -- 본 분석은 pk_parameters_v3.R의 0.35 저자가정 사용)",
    model       = "1-compartment",
    calibration_note = "Paxil 2017 label 30 mg/day 30일 정상상태 Cmax 61.7, Cmin 30.7 ng/mL, Tmax 5.2h, t_half 21h에 보정",
    source      = "Paxil 2017 FDA label (Cmax/Cmin/Tmax/t_half 보정); Li 2022 [8] Table 2 (CV_CL -- 중국 환자 184명); CV_Vd는 저자 가정 35%(추가 지시 1 2-1, Li 2022 V/F는 misidentified로 판단해 미사용)"
  ),

  # -- Venlafaxine IR -------------------------------------------
  # D7: parent, ODV CL/F를 XR 150mg QD 라벨 정상상태 Cmax(150/260 ng/mL)에
  # brentq로 보정 (calibration/ven_calibration_v3.py). V/F는 라벨값 고정.
  # IR의 Ka만 별도로 라벨 Tmax에서 유도 -- IR 투여설계(25mg TID, 본 분석)는
  # 독립 점검용이며, parent/ODV CL/F, V/F 자체는 XR과 공유 (같은 약물,
  # 같은 전신 청소율을 가정).
  venlafaxine_IR = list(
    drug        = "Venlafaxine",
    class       = "SNRI",
    formulation = "IR",
    Ka          = 1.072,     # 1/h | IR 75mg BID 라벨 정상상태 Tmax 2h에서 유도
    Ka_basis    = "label_derived",
    Vd          = 525,       # L (V/F) | 7.5 L/kg x 70kg | Effexor XR 2017 label (고정)
    Vd_basis    = "label",
    CL          = 66.3,      # L/h (CL/F) = 0.95 L/h/kg | D7: XR 150mg QD 라벨
                              # Cmax 150에 brentq 보정
    CL_basis    = "label_calibrated (D7, brentq to XR label Cmax)",
    F           = 1.0,
    F_basis     = "apparent (CL/F, V/F로 해석)",
    t_half      = 5.5,       # h | V/F, CL/F에서 유도
    CV_CL       = 0.46,      # 라벨 SD/평균 1.3+-0.6 L/h/kg (비율은 재보정
                              # 전후 동일하게 적용)
    CV_Vd       = 0.49,      # 라벨 SD/평균 7.5+-3.7 L/kg
    model       = "1-compartment + metabolite (ODV)",
    ODV = list(
      Vd_odv       = 399,    # L | 5.7 L/kg x 70kg | Effexor XR 2017 label (고정)
      Vd_odv_basis = "label",
      CL_odv       = 21.1,   # L/h = 0.30 L/h/kg | D7: XR 라벨 ODV Cmax 260에
                              # brentq 보정, Fm=0.80 조건부
      CL_odv_basis = "label_calibrated (D7, brentq to XR label ODV Cmax; conditional on Fm=0.80)",
      t_half       = 13.1,   # h
      CV_CL_odv    = 0.50,   # 라벨 SD/평균 0.4+-0.2 L/h/kg
      Fm           = 0.80,
      Fm_basis     = "unsourced_pending (author assumption) -- ODV CL/F 보정은 이 Fm 값 조건부"
    ),
    source = "Effexor XR 2017 FDA label (V/F 고정); calibration/ven_calibration_v3.py (CL/F brentq 보정, D7)",
    note   = "venlafaxine IR 투여설계(25 mg TID)는 독립 점검용 -- 본 분석의 CL/F, V/F 보정은 XR 라벨 기준"
  ),

  # -- Venlafaxine XR -------------------------------------------
  venlafaxine_XR = list(
    drug        = "Venlafaxine",
    class       = "SNRI",
    formulation = "XR",
    Ka          = 0.216,     # 1/h | XR 150mg QD 라벨 정상상태 Tmax 5.5h에서 유도
    Ka_basis    = "label_derived",
    Vd          = 525,
    Vd_basis    = "label",
    CL          = 66.3,
    CL_basis    = "label_calibrated (D7, brentq to XR label Cmax)",
    F           = 1.0,
    F_basis     = "apparent (CL/F, V/F로 해석)",
    t_half      = 5.5,
    CV_CL       = 0.46,
    CV_Vd       = 0.49,
    model       = "1-compartment + metabolite (ODV)",
    ODV = list(
      Vd_odv       = 399,
      Vd_odv_basis = "label",
      CL_odv       = 21.1,
      CL_odv_basis = "label_calibrated (D7, brentq to XR label ODV Cmax; conditional on Fm=0.80)",
      t_half       = 13.1,
      CV_CL_odv    = 0.50,
      Fm           = 0.80,
      Fm_basis     = "unsourced_pending (author assumption) -- ODV CL/F 보정은 이 Fm 값 조건부"
    ),
    source = "Effexor XR 2017 FDA label (V/F, Cmax/ODV Cmax 보정목표); calibration/ven_calibration_v3.py (CL/F brentq 보정, D7)"
  )
)

# -- CYP2D6 genotype 보정 (D5) --------------------------------
# v3: paroxetine에만 적용. PM x 0.58 (Sindrup 1992, Clin Pharmacol Ther
# 51:278-87, Table II -- 정상상태 30mg/day 총경구 청소율 중앙값
# EM 36 L/h / PM 21 L/h = 0.58; AUC비 1.73의 역수와 동일).
# RM 범주는 삭제 (근거수준 부족, evidence/CYP2D6_EVIDENCE_1027.md 3절:
# "RM(1.60, 빈도22%) ... 빼는 것을 권합니다"). fluoxetine, venlafaxine에는
# phenotype을 배정하지 않음 (D5, 모델링 안 함).
#
# PM 빈도 8%: Gaedigk, Sangkuhl, Whirl-Carrillo et al. "Prediction of
# CYP2D6 phenotype from genotype across world populations." Genet Med.
# 2017;19:69-76. Discussion 문단의 유럽인 "관찰된 표현형(phenotypic)
# PM 빈도(pPM)" 8.45%와 일치 (2026-10-28 WebFetch로 재확인). Table 1/
# Supplementary Table S2의 유전형 예측치(유럽인 평균 5.4%)와는 다른
# 수치이므로, 원고에는 "phenotypic PM frequency (pPM), Discussion"로
# 명시할 것 -- Table 1 (genotype-predicted)을 인용하지 말 것.
#
# normal_metabolizer freq = 1 - PM freq: RM 범주를 삭제하면서 paroxetine
# 집단이 PM/비PM 두 범주만 남음에 따른 필연적 귀결 (HANDOFF_v3가 직접
# 지정한 숫자는 아니지만, 7-2 기대값 섹션의 "PM 8%일 때 (NM 중앙값 대비)
# 약 4% 길어짐"이라는 서술과 정확히 부합함 -- PM 8%/비PM 92%로 가중평균한
# 결과와 일치하는지 7-2에서 재확인할 것).
cyp2d6 <- list(
  poor_metabolizer = list(
    freq            = 0.08,
    CL_factor       = 0.58,
    freq_basis      = "Gaedigk 2017 [15] Discussion -- phenotypic PM frequency (pPM), Europeans, 8.45% (NOT Table 1 genotype-predicted 5.4%)",
    CL_factor_basis = "Sindrup 1992, Clin Pharmacol Ther 51:278-87, Table II -- steady-state 30 mg/day total oral clearance ratio PM/EM = 21/36 = 0.58"
  ),
  normal_metabolizer = list(
    freq            = 0.92,
    CL_factor       = 1.00,
    freq_basis      = "1 - PM frequency (RM 범주 삭제로 단일 reference 범주가 됨, D5)",
    CL_factor_basis = "reference (배수 1.0)"
  )
)
# rapid_metabolizer (RM) 범주는 v3에서 완전히 삭제됨 -- HANDOFF_v3 5절
# "삭제할 것: ... CYP2D6 RM 범주"; evidence/CYP2D6_EVIDENCE_1027.md 3절:
# "RM(1.60, 빈도 22%): CPIC 범주는 PM/IM/NM/UM이고 ... 22%라는 빈도도
# UM 빈도로 보기 어렵습니다. 빼는 것을 권합니다."

# sensitivity_v3.R 전용 -- 주 분석에는 쓰지 않음 (HANDOFF_v3 7-1-4).
# paroxetine PM 배수 민감도: 0.58 (주분석, Sindrup) vs 0.40 (Liao 2024
# 보정 Css비 2.50의 역수).
cyp2d6_pm_sensitivity_alt <- list(
  CL_factor       = 0.40,
  CL_factor_basis = "Liao 2024 -- calibrated Css ratio 2.50의 역수 (sensitivity_v3.R 전용, 주분석 아님)"
)

# -- 검증 요약 --------------------------------------------------
cat("*** parV835 변형 로드됨 -- paroxetine CV_Vd=0.835 (게이트 대조 전용, 본 분석 아님) ***\n")
cat("PK 파라미터 v3 로드 완료 (HANDOFF_TO_NARI_v3.md 5절 Table 1 반영)\n")
cat("약물:", paste(names(pk_params), collapse=" | "), "\n\n")
cat("[v3 변경사항 -- v2 대비]\n")
cat("  sertraline: Vd 1400->3054 L, CL 37.3->77.8 L/h (D2', Zoloft 2009 label AUC 재보정)\n")
cat("  paroxetine: Vd 609->813 L, CL 20.1->26.8 L/h, Ka 0.580->0.39 (D2', Paxil 2017 label Cmax/Cmin 재보정)\n")
cat("  paroxetine CV: CL 60%->51.9% (D9, Li 2022 Table 2); Vd 83.5% (parV835 -- 본 분석은 0.35, pk_parameters_v3.R 참고)\n")
cat("  venlafaxine: CL 91->66.3 L/h, ODV CL 28->21.1 L/h, t_half 5->5.5h / 11->13.1h (D7, XR label Cmax brentq 보정)\n")
cat("  fluoxetine: 값 유지, 출처만 NZ Medicines -> Christchurch Medicines Information Service (1999)\n")
cat("  norfluoxetine Ka_met 삭제 (모델식에 대사체 Ka 없음, B4)\n")
cat("  CYP2D6: paroxetine에만 적용 (PM x0.58, Sindrup 1992), RM 범주 삭제 (D5)\n")
cat("  CYP2D6 PM 빈도 8% 출처 재확인: Gaedigk 2017 Discussion (phenotypic, 8.45%) -- Table 1 (genotype, 5.4%) 아님\n")
