%% 1. 초기 상수 정의
% clc; clear all; close all; format longG;

% 지구 상수 (SI 단위, m/kg/s)
mu = 3.986004418e14; % (m^3/s^2)
R_e = 6378.137e3;    % 지구 반경 (m)
J2 = 0.00108263;     % J2 계수

%% 2. 초기 궤도 요소 (과제 조건)
a0_km = 6878;
e0 = 0.01;
i0_deg = 40;
RAAN0_deg = 120;
AoP0_deg = 20;
M0_deg = 60;

% --- 계산을 위해 단위를 (m, rad)로 변환 ---
a0 = a0_km * 1000;         % m
i0 = deg2rad(i0_deg);      % rad
RAAN0 = deg2rad(RAAN0_deg);% rad
AoP0 = deg2rad(AoP0_deg);  % rad
M0 = deg2rad(M0_deg);      % rad

%% 3. 영년 변화율(Secular Rates) 및 시간 벡터
% --- 중간 변수 계산 ---
p0 = a0 * (1 - e0^2);         % Semi-latus rectum (m)
n0 = sqrt(mu / a0^3);         % 평균 운동 (rad/s) 의 n

% --- 궤도 요소별 영년 변화율 (rad/s) ---
a_rate = 0; 
e_rate = 0; 
i_rate = 0; 
RAAN_rate = -(3/2) * n0 * J2 * (R_e / p0)^2 * cos(i0); 
AoP_rate = (3/4) * n0 * J2 * (R_e / p0)^2 * (5 * cos(i0)^2 - 1); 
M_perturb_rate = (3/4) * n0 * J2 * (R_e / p0)^2 * sqrt(1 - e0^2) * (3 * cos(i0)^2 - 1); 
M_rate_total = n0 + M_perturb_rate;

% --- 시간 벡터 생성 (1일) ---
t_end_days = 1/8;
t_end_sec = t_end_days * 86400;
dt_sec = 60; % 1분 간격
t_vec = (0:dt_sec:t_end_sec)';
num_steps = length(t_vec);
t_vec_days = t_vec / 86400;

%% 4. "평균" 궤도 요소 전파 (Secular Mean Elements)
a_hist = ones(num_steps, 1) * a0;
e_hist = ones(num_steps, 1) * e0;
i_hist = ones(num_steps, 1) * i0;
p_hist = ones(num_steps, 1) * p0;
RAAN_hist = RAAN0 + RAAN_rate * t_vec;
AoP_hist = AoP0 + AoP_rate * t_vec;
M_hist = M0 + M_rate_total * t_vec;

%% 5. 단주기 계산을 위한 입력 변수 계산
nu_hist = zeros(num_steps, 1);
r_hist = zeros(num_steps, 1);

fprintf('Calculating mean nu and r histories...\n');
for k = 1:num_steps
    nu_scalar = M2nu(M_hist(k), e_hist(k),i_hist(k)); 
    nu_hist(k) = nu_scalar;
    r_hist(k) = p_hist(k) / (1 + e_hist(k) * cos(nu_scalar));
end

%% 6. 헬퍼 변수 계산 (Loop 이전)
j = 2;
temp_e_term = (-e0) / (1 + sqrt(1 - e0^2));
cos2v_bar = (temp_e_term^j) * (1 + j * sqrt(1 - e0^2));

%% 7. 최종 순간 궤도 요소(Osculating Elements) 계산
a_osc_hist = zeros(num_steps, 1);
e_osc_hist = zeros(num_steps, 1);
i_osc_hist = zeros(num_steps, 1);
RAAN_osc_hist = zeros(num_steps, 1);
AoP_osc_hist = zeros(num_steps, 1);
M_osc_hist = zeros(num_steps, 1);
M0_osc_hist = zeros(num_steps,1);

fprintf('Calculating full short-periodic corrections...\n');

for k = 1:num_steps
    % --- k번째 시점의 '평균' 값 가져오기 ---
    a_k     = a_hist(k);
    e_k     = e_hist(k);
    i_k     = i_hist(k);
    p_k     = p_hist(k);
    r_k     = r_hist(k);
    omega_k = AoP_hist(k);
    nu_k    = nu_hist(k);
    M_k     = M_hist(k); % M_k는 n0*t를 포함한 전체 평균 근점 이각입니다.
    
    % --- (A) 1주기 '총 변화량'(Delta) 계산 ---
    
    delta_a_total = (J2 * R_e^2 / a_k) * ...
        ( (a_k/r_k)^3 - 1/(1 - e_k^2)^(3/2) + ...
        ( -(a_k/r_k)^3 + 1/(1 - e_k^2)^(3/2) + (a_k/r_k)^3 * cos(2*omega_k + 2*nu_k) ) ...
        * (3 * sin(i_k)^2 / 2) );
    
    delta_e_total = (J2 * R_e^2 / 4) * ...
        ( ...
        (-2) / (a_k^2 * e_k * sqrt(1 - e_k^2)) ...
        + 2 * a_k * (1 - e_k^2) / (e_k * r_k^3) ...
        + ( ...
            3 / (a_k^2 * e_k * sqrt(1 - e_k^2)) ...
            - 3 * a_k * (1 - e_k^2) / (e_k * r_k^3) * cos(2*omega_k + nu_k) ...
            - 3 * (1 - e_k^2) / p_k^2 * cos(2*omega_k + 2*nu_k) ...
            - 3 * cos(2*omega_k + 2*nu_k) / (a_k^2 * e_k * (1 - e_k^2)) ...
            + 3 * a_k * (1 - e_k^2) / (e_k * r_k^3) * cos(2*omega_k + 2*nu_k) ...
            - (1 - e_k^2) / p_k^2 * cos(2*omega_k + 3*nu_k) ...
          ) ...
        ) * sin(i_k)^2;

    delta_inc_total = (J2 * R_e^2 * sin(2*i_k)) / (8 * p_k^2) * ...
        ( 3 * cos(2*omega_k + 2*nu_k) + 3 * e_k * cos(2*omega_k + nu_k) + e_k * cos(2*omega_k + 3*nu_k) );

    E_k = M2E(M_k, e_k);
    v_M_diff = E_k - M_k; 
    
    delta_RAAN_total = -(J2 * R_e^2 * cos(i_k)) / (4 * p_k^2) * ...
        ( 6 * (v_M_diff) - 3 * sin(2*omega_k + 2*nu_k) - 3 * e_k * sin(2*omega_k + nu_k) - e_k * sin(2*omega_k + 3*nu_k) );

    delta_omega_total = (3 * J2 * R_e^2) / (2 * p_k^2) * ( ...
          (2 - (5/2) * sin(i_k)^2) * (v_M_diff) ... % (nu - M + e*sin(nu)) 대체
        + (1 - (3/2) * sin(i_k)^2) * ( (1/e_k) * (1 - 0.25 * e_k^2) * sin(nu_k) + 0.5 * sin(2*nu_k) + (e_k/12) * sin(3*nu_k) ) ...
        - (1/e_k) * ( 0.25 * sin(i_k)^2 + (0.5 - 15/16 * sin(i_k)^2) * e_k^2 ) * sin(2*omega_k + nu_k) ...
        + (e_k/16) * sin(i_k)^2 * sin(nu_k - 2*omega_k) - 0.5 * (1 - 2.5 * sin(i_k)^2) * sin(2*omega_k + 2*nu_k) ...
        + (1/e_k) * ( (7/12) * sin(i_k)^2 - (1/6) * (1 - 19/8 * sin(i_k)^2) * e_k^2 ) * sin(2*omega_k + 3*nu_k) ...
        + (3/8) * sin(i_k)^2 * sin(4*nu_k + 2*omega_k) + (e_k/16) * sin(i_k)^2 * sin(2*omega_k + 5*nu_k) );

    delta_M_total = (3 * J2 * R_e^2 * sqrt(1 - e_k^2)) / (2 * e_k * p_k^2) * ( ...
        -(1 - 1.5 * sin(i_k)^2) * ( (1 - 0.25 * e_k^2) * sin(nu_k) + (e_k/2) * sin(2*nu_k) + (e_k^2/12) * sin(3*nu_k) ) ...
        + sin(i_k)^2 * ( 0.25 * (1 + 5/4 * e_k^2) * sin(2*omega_k + nu_k) - (e_k^2/16) * sin(nu_k - 2*omega_k) - (7/12) * (1 - e_k^2/28) * sin(2*omega_k + 3*nu_k) ) ...
        - (3 * e_k / 8) * sin(2*omega_k + 4*nu_k) - (e_k^2 / 16) * sin(2*omega_k + 5*nu_k) );

    
    % --- (B) 1주기 '평균값'(Delta_bar) 계산 ---
    
    delta_a_mean = 0; 

    mean_e_factor = (1/4) * (J2 * R_e^2 / p_k^2) * sin(i_k)^2 * ((1-e_k^2)/e_k) * cos(2*omega_k);
    delta_e_mean = mean_e_factor * cos2v_bar; 
    
    mean_inc_factor = -(1/8) * (J2 * R_e^2 / p_k^2) * sin(2*i_k)^2 * cos(2*omega_k);
    delta_inc_mean = mean_inc_factor * cos2v_bar; 

    mean_RAAN_factor = -(1/4) * (J2 * R_e^2 / p_k^2) * cos(i_k) * sin(2*omega_k);
    delta_RAAN_mean = mean_RAAN_factor * cos2v_bar; 

    term1_omega = sin(i_k)^2 * (1/8 + (1-e_k^2)/(6*e_k^2) * cos2v_bar);
    term2_omega = (1/6) * cos(i_k)^2 * cos2v_bar;
    delta_omega_mean = (3/2) * (J2 * R_e^2 / p_k^2) * (term1_omega + term2_omega) * sin(2*omega_k); 

    term1_M = (1/8) + (1 + e_k^2/2) / (6*e_k^2) * cos2v_bar;
    delta_M_mean = -(3/2 * J2 * R_e^2 / p_k^2) * sqrt(1 - e_k^2) * sin(i_k)^2 * term1_M * sin(2*omega_k); 

    
    % --- (C) 최종 순간 궤도 요소(Osculating) 계산 ---
    a_osc_hist(k)    = a_hist(k)    + (delta_a_total    - delta_a_mean);
    e_osc_hist(k)    = e_hist(k)    + (delta_e_total    - delta_e_mean);
    i_osc_hist(k)    = i_hist(k)    + (delta_inc_total  - delta_inc_mean);
    RAAN_osc_hist(k) = RAAN_hist(k) + (delta_RAAN_total - delta_RAAN_mean);
    AoP_osc_hist(k)  = AoP_hist(k)  + (delta_omega_total - delta_omega_mean);
    M_osc_hist(k)    = M_hist(k)    + (delta_M_total - delta_M_mean); 
    
end
fprintf('Calculation complete.\n');

%% 8. 최종 결과 시각화 (Osculating Elements)
figure('Name', 'Analytical Short-Periodic Perturbations (J2)');
sgtitle('Osculating Elements (Secular + Short-Periodic)', 'FontSize', 14, 'FontWeight', 'bold');

subplot(3, 2, 1);
plot(t_vec_days, a_osc_hist / 1000); % m -> km
title('Semi-major Axis (a)'); xlabel('Time (days)'); ylabel('km'); grid on;

subplot(3, 2, 2);
plot(t_vec_days, e_osc_hist);
title('Eccentricity (e)'); xlabel('Time (days)'); ylabel('Unitless'); grid on;

subplot(3, 2, 3);
plot(t_vec_days, rad2deg(i_osc_hist)); % rad -> deg
title('Inclination (i)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;

subplot(3, 2, 4);
plot(t_vec_days, rad2deg(mod(RAAN_osc_hist, 2*pi))); % Wrap 0-360
title('RAAN (\Omega) (wrapped)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;

subplot(3, 2, 5);
plot(t_vec_days, rad2deg(mod(AoP_osc_hist, 2*pi))); % Wrap 0-360
title('Argument of Perigee (\omega) (wrapped)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;

subplot(3, 2, 6);
plot(t_vec_days, rad2deg(mod(M_osc_hist, 2*pi))); % Wrap 0-360
title('Mean Anomaly (M) (wrapped)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;

%% 9. [신규] 영년 효과가 제거된 순수 단주기 효과 시각화
% figure('Name', 'Pure Short-Periodic Oscillations (Delta - Delta_bar)');
% sgtitle('Osculating Elements minus Secular Trend (Pure Periodic Variations)', 'FontSize', 14, 'FontWeight', 'bold');

% 계산: Osculating History - Mean History
a_pure_sp = a_osc_hist - a_hist;
e_pure_sp = e_osc_hist - e_hist;
i_pure_sp = i_osc_hist - i_hist;
RAAN_pure_sp = RAAN_osc_hist - RAAN_hist;
AoP_pure_sp = AoP_osc_hist - AoP_hist;
M_pure_sp = M_osc_hist - M_hist; % M의 섭동 부분만 (n*t 제외)

subplot(3, 2, 1);
plot(t_vec_days, a_pure_sp / 1000); % m -> km
title('Detrended Semi-major Axis (a)'); xlabel('Time (days)'); ylabel('km'); grid on;

subplot(3, 2, 2);
plot(t_vec_days, e_pure_sp);
title('Detrended Eccentricity (e)'); xlabel('Time (days)'); ylabel('Unitless'); grid on;

subplot(3, 2, 3);
plot(t_vec_days, rad2deg(i_pure_sp)); % rad -> deg
title('Detrended Inclination (i)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;

subplot(3, 2, 4);
plot(t_vec_days, rad2deg(RAAN_pure_sp)); % rad -> deg
title('Detrended RAAN (\Omega)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;

subplot(3, 2, 5);
plot(t_vec_days, rad2deg(AoP_pure_sp)); % rad -> deg
title('Detrended Argument of Perigee (\omega)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;

subplot(3, 2, 6);
plot(t_vec_days, rad2deg(M_pure_sp)); % rad -> deg
title('Detrended Mean Anomaly (M0 correction)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;

function E = M2E(M, e)

    options = optimset('Display','off', 'TolFun', 1e-12);
    % 케플러 방정식: E - e*sin(E) - M = 0
    E = fsolve(@(E_guess) E_guess - e*sin(E_guess) - M, M, options);
end