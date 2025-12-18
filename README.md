# J2 Perturbation Study (Analytical vs Numerical)

이 레포지토리는 지구 비구면성(특히 **J2**)에 의한 궤도 요소 변화(영년/단주기)를  
**해석적(Analytical) 모델**과 **수치적(Cowell) 전파**로 비교·분석하기 위한 MATLAB 코드 모음입니다.

- 대상 궤도(과제 조건 예시):  
  a = 6878 km, e = 0.01 (또는 e = 0), i = 40°, RAAN = 120°, AoP = 20°, M0 = 60°
- 주요 산출물:  
  - 영년 변화율(Ω̇, ω̇, Ṁ 등)  
  - 단주기(Short-period) 섭동이 포함된 순간(oscillating) 궤도요소  
  - 수치 전파 결과의 추세(Trend) 추정 및 detrend(진동 분리) 분석

---

## 파일 구성 및 역할

### 1) `MAIN_J2_Analytical.m`
**J2 영년(Secular) 변화율**만을 사용해 궤도요소를 해석적으로 전파합니다.

- 하는 일
  - 지구 상수(μ, Re, J2) 및 초기 궤도요소 입력
  - J2 기반 영년 변화율 계산:
    - ȧ = 0, ė = 0, i̇ = 0
    - Ω̇(승교점 적경), ω̇(근지점 인수), Ṁ(평균근점이각)의 영년 변화율
  - `Element(t) = Element0 + Rate * t` 형태로 시간 전파
  - 1분 간격으로 시간 히스토리 생성 및 플롯
  - 1일 후 변화량/변화율을 로그로 출력

- 출력(시각화)
  - RAAN(Ω), AoP(ω), Mean Anomaly(M), M0(섭동항만 분리한 값) 시간 변화 그래프

> 참고: 본 스크립트는 “영년항”만 포함하므로, 단주기 진동은 나타나지 않습니다.

---

### 2) `MAIN_J2_Analytical_ShortPeriod.m`
J2의 **영년 + 단주기(Short-period) 보정**을 포함한 해석적 “순간(oscillating)” 궤도요소를 계산합니다.

- 하는 일(흐름)
  1. 영년 변화율을 이용해 **평균(mean) 궤도요소** 전파  
     (`RAAN_hist`, `AoP_hist`, `M_hist` 등)
  2. 각 시간스텝에서 평균근점이각 M → 편심이각 E → 진근점이각 ν(또는 ν 계산 함수)로 변환하고 r 계산
  3. Lyddane/고전 해석식 형태의 단주기 보정량 Δ(총 변화) 및 평균 성분 Δ̄(평균값) 계산
  4. **Osculating = Mean + (Δ − Δ̄)** 로 순간 궤도요소 히스토리 생성
  5. Osculating 요소 플롯 + “추세 제거(osc-mean)”에 해당하는 순수 단주기 진동도 플롯

- 출력(시각화)
  - (Secular+Short-period) Osculating 요소: a, e, i, Ω, ω, M (wrapped)
  - “순수 단주기 성분”: (Osculating − Mean) 형태의 detrended 진동

- 의존/헬퍼 함수
  - `M2nu(...)` : M → ν 변환(외부 함수로 존재해야 함)
  - `M2E(...)` : 내부에 정의(Kepler 방정식 fsolve 사용)

> 주의: 본 스크립트는 단주기 항이 많아 수식/부호/각도 wrap 여부(unwrap/mod)에 민감합니다.  
> 결과 해석 시 Ω, ω, M의 wrap 처리 기준을 일관되게 유지하세요.

---

### 3) `MAIN_J2_Analytical_ShortPeriod_Circular.m`
원궤도(e = 0) 가정 하에서, **J2 단주기 섭동(원궤도 근사식)**을 사용해 Osculating 요소를 계산합니다.

- 특징
  - e=0이므로 ω와 M 대신 **위도인수 u(= ω + ν)**를 사용
  - 영년 변화:
    - Ω̇ 계산
    - u̇ = n + (J2로 인한 Ṁ 섭동항 + ω̇) 형태로 근사 구성
  - 단주기 보정:
    - Δa, Δi, ΔΩ, Δu를 u(평균 위도인수)의 sin(2u), cos(2u) 형태로 계산
  - Osculating = Mean + Δ

- 출력(시각화)
  - Osculating: a, i, Ω, u (wrap)
  - 부록 플롯: (Osculating − Mean) = 순수 단주기 진동(Δ) 확인

> 이 파일은 “원궤도 근사”가 핵심이므로, e≈0 조건에서 단주기 특성만 빠르게 확인하는 용도로 적합합니다.

---

### 4) `MAIN_J2_Numerical.m`
J2만 포함한 동역학 모델로 **수치적(Cowell) 전파**를 수행하고,  
전파 결과로부터 **순간 궤도요소(Osculating)**를 추출한 뒤 **영년 추세/진동 성분**을 분리합니다.

- 하는 일(핵심)
  - 전파 환경 설정:
    - `global_variable`, `STEP1_constants()` 등 사용자 정의 환경 사용
    - 섭동 모델: J2만 활성화(in.n=2, 나머지(태양/달/drag 등) off)
  - 초기 궤도요소 → 초기 상태벡터 변환 후 ODE 적분(`ode89`, `Accel`)
  - 상태벡터 히스토리 → 궤도요소 히스토리 변환(`processOrbitData`)
  - 영년 추세(Secular trend) 추정:
    - 각 궤도요소에 대해 polyfit(선형)으로 trend 산출(각도는 unwrap 후 fitting)
  - “진동 성분” 분리:
    - detrended = osculating − trend
  - 30일 기준으로 평균 변화율(deg/day)을 계산하고 출력

- 출력(시각화)
  1) Osculating 요소 시간변화 플롯 (a, e, i, Ω, ω, ν, M, M0 등)  
  2) Osculating vs Trend 비교 플롯  
  3) detrended(osc − trend) 플롯: 순수 주기 성분 확인  
  4) 30일 기준 영년 변화율 로그 출력

- 의존 함수/파일(레포에 포함되어 있어야 함)
  - `global_variable`, `STEP1_constants`
  - `Mjday`, `Accel`, `ode89`
  - `coe2rv`, `M2nu`, `processOrbitData`

> 주의: 이 스크립트는 사용자 프레임워크(상수/가속도/ODE/요소변환)에 강하게 의존합니다.  
> 레포 공개 시에는 해당 의존 파일을 함께 포함하거나, 사용 방법을 별도로 명시하는 것을 권장합니다.

---

### 5) `MAIN_twobody.m`
가장 단순한 **2체(central gravity only)** 전파를 수행해,  
수치 적분의 기본 성질(특히 에너지 보존)을 확인하는 검증용 스크립트입니다.

- 하는 일
  - 2체 ODE: ṙ = v, v̇ = −μ r/|r|^3
  - `ode89`로 1일 전파
  - 특정 기계적 에너지 ξ = v^2/2 − μ/r 계산 후 시간에 따른 변화(보존성) 확인
  - (선택) 상태벡터를 궤도요소로 변환하여 요소 변화 플롯(`rv2coe`)

- 출력(시각화)
  - Specific mechanical energy(ξ) vs time
  - (선택) a, e, i, Ω, ω, ν/M 등의 요소 플롯

> 참고: 2체 전파에서는 (수치 오차를 제외하면) 에너지가 일정해야 합니다.  
> 본 파일은 이후 J2 수치전파의 결과 신뢰성을 점검하는 baseline으로 활용할 수 있습니다.

---

## 실행 가이드(권장 순서)

1. **해석적 영년만 확인**: `MAIN_J2_Analytical.m`  
2. **해석적 단주기 포함**: `MAIN_J2_Analytical_ShortPeriod.m`  
3. **원궤도 단주기 근사**: `MAIN_J2_Analytical_ShortPeriod_Circular.m`  
4. **수치 전파 기반 비교**: `MAIN_J2_Numerical.m`  
5. **2체 baseline 검증**: `MAIN_twobody.m`

---

## 비고 / 주의사항

- 각도 요소(Ω, ω, M, ν)는 분석/피팅 시 **unwrap 처리** 후 변화율을 구하고,  
  플롯에서는 필요 시 `mod(·, 2π)`로 wrap하여 표시하는 것이 혼동을 줄입니다.
- `MAIN_J2_Numerical.m`은 외부 의존성이 많으므로, 레포에 포함된 함수 목록을 확인하세요.
****
