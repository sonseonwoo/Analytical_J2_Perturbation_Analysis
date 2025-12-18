%% 1. 초기 상수 정의
clc; clear all; close all; format longG;

% 지구 상수 (WGS84 기준, SI 단위)
mu = 3.986004418e14; % (m^3/s^2)
R_e = 6378.137e3;    % 지구 반경 (m)
J2 = 0.00108263;     % J2 계수

%% 2. 초기 궤도 요소 (원궤도 가정)
% --- 과제 조건 (이심률 e=0 설정) ---
a0_km = 6878;
e0 = 0;              % [중요] 원궤도 가정
i0_deg = 40;
RAAN0_deg = 120;
u0_deg = 80;         % 원궤도에서는 M과 w 대신 위도인수(u = w+M)를 사용

% --- 계산을 위해 단위를 (m, rad)로 변환 ---
a0 = a0_km * 1000;         % m
i0 = deg2rad(i0_deg);      % rad
RAAN0 = deg2rad(RAAN0_deg);% rad
u0 = deg2rad(u0_deg);      % rad (Argument of Latitude)

%% 3. 영년 변화율(Secular Rates) 및 시간 벡터
% --- 중간 변수 계산 ---
% 원궤도이므로 p = a, n = sqrt(mu/a^3)
p0 = a0; 
n0 = sqrt(mu / a0^3);         % 평균 운동 (rad/s)
p_factor = (R_e / p0)^2;

% --- 궤도 요소별 영년 변화율 (rad/s) ---
% a, e, i는 J2 섭동에 의한 영년 변화가 없음 (0)
a_rate = 0;
i_rate = 0;

% RAAN(승교점 적경)의 영년 변화율
RAAN_rate = -1.5 * n0 * J2 * p_factor * cos(i0);

% u(위도 인수)의 영년 변화율 (u_dot = n + M_dot + w_dot)
% e=0 일 때의 근사식 적용
M_w_rate_term = 0.75 * n0 * J2 * p_factor;
w_dot = M_w_rate_term * (4 - 5*sin(i0)^2);
M_dot_perturb = M_w_rate_term * (2 - 3*sin(i0)^2); % n0 제외한 섭동항만
u_rate = n0 + M_dot_perturb + w_dot; 

% --- 시간 벡터 생성 (1일, 1분 간격) ---
t_end_days = 1.0;          % 1일
t_end_sec = t_end_days * 86400;
dt_sec = 60;               % 60초 간격
t_vec = (0:dt_sec:t_end_sec)';
num_steps = length(t_vec);
t_vec_days = t_vec / 86400;

%% 4. "평균" 궤도 요소 전파 (Secular Mean Elements)
% a, i, e는 상수, 각도(RAAN, u)는 선형 증가
a_mean_hist = ones(num_steps, 1) * a0;
i_mean_hist = ones(num_steps, 1) * i0;
RAAN_mean_hist = RAAN0 + RAAN_rate * t_vec;
u_mean_hist = u0 + u_rate * t_vec;

%% 5. 단주기 섭동 계산 (Short-Periodic Perturbations - Lyddane/SGP4)
% 최종 순간 궤도 요소(Osculating Elements) 저장용 배열
a_osc_hist = zeros(num_steps, 1);
i_osc_hist = zeros(num_steps, 1);
RAAN_osc_hist = zeros(num_steps, 1);
u_osc_hist = zeros(num_steps, 1);

fprintf('Calculating full short-periodic corrections (Circular Case)...\n');

for k = 1:num_steps
    % --- k번째 시점의 '평균' 값 가져오기 ---
    a_k    = a_mean_hist(k);
    i_k    = i_mean_hist(k);
    u_k    = u_mean_hist(k); % 평균 위도 인수
    raan_k = RAAN_mean_hist(k);
    
    % --- 공통 계산 항 (Lyddane 1963 / SGP4 Circular) ---
    sin_i = sin(i_k);
    cos_i = cos(i_k);
    sin2u = sin(2 * u_k);
    cos2u = cos(2 * u_k);
    
    % J2 계수 항
    temp1 = 1.5 * J2 * (R_e / a_k)^2;
    temp2 = 0.5 * temp1; % = 0.75 * J2 * ...
    
    % --- (A) 단주기 보정량(Delta) 계산 ---
    % 1. 장반경 (Semi-major Axis) 섭동
    % 식: (3*J2*Re^2 / 4a) * [sin^2(i)*cos(2u) + (3cos^2(i)-1)]
    % 코드의 temp1*a_k/2 = (3/2 * J2 * Re^2/a) / 2 = 3/4 * ... (계수 일치)
    term_a = (sin_i^2) * cos2u + (3 * cos_i^2 - 1);
    delta_a = (a_k * temp1 / 2) * term_a; 
    
    % 2. 경사각 (Inclination) 섭동
    % 식: (9*J2*Re^2 / 8a^2) * cos(i)*sin(i)*cos(2u)
    delta_i = 1.5 * temp2 * cos_i * sin_i * cos2u;
    
    % 3. 승교점 적경 (RAAN) 섭동
    % 식: (9*J2*Re^2 / 8a^2) * cos(i)*sin(2u)
    delta_raan = 1.5 * temp2 * cos_i * sin2u;
    
    % 4. 위도 인수 (Argument of Latitude) 섭동
    % 식: -(3*J2*Re^2 / 16a^2) * (7cos^2(i)-1) * sin(2u)
    delta_u = -0.25 * temp2 * (7 * cos_i^2 - 1) * sin2u;
    
    % --- (B) 최종 순간 궤도 요소(Osculating) 계산 ---
    a_osc_hist(k)    = a_k + delta_a;
    i_osc_hist(k)    = i_k + delta_i;
    RAAN_osc_hist(k) = raan_k + delta_raan;
    u_osc_hist(k)    = u_k + delta_u;
    
end
fprintf('Calculation complete.\n');

%% 6. 최종 결과 시각화 (Osculating Elements)
figure('Name', 'Analytical Short-Periodic Perturbations (Circular J2)', 'Color', 'w');
sgtitle('Osculating Elements (Circular Orbit Approximation)', 'FontSize', 14, 'FontWeight', 'bold');

% 1) Semi-major Axis
subplot(2, 2, 1);
plot(t_vec_days, a_osc_hist / 1000, 'LineWidth', 1.2); % m -> km
title('Semi-major Axis (a)'); 
xlabel('Time (days)'); ylabel('km'); grid on; xlim([0 t_end_days]);

% 2) Inclination
subplot(2, 2, 2);
plot(t_vec_days, rad2deg(i_osc_hist), 'LineWidth', 1.2); % rad -> deg
title('Inclination (i)'); 
xlabel('Time (days)'); ylabel('Degrees'); grid on; xlim([0 t_end_days]);

% 3) RAAN (Wrapped)
subplot(2, 2, 3);
plot(t_vec_days, rad2deg(mod(RAAN_osc_hist, 2*pi)), 'LineWidth', 1.2);
title('RAAN (\Omega)'); 
xlabel('Time (days)'); ylabel('Degrees'); grid on; xlim([0 t_end_days]);

% 4) Argument of Latitude (Wrapped)
subplot(2, 2, 4);
% u는 빠르게 증가하므로 360도로 mod 연산하여 표현
plot(t_vec_days, rad2deg(mod(u_osc_hist, 2*pi)), 'LineWidth', 1.2);
title('Argument of Latitude (u = \omega + \nu)'); 
xlabel('Time (days)'); ylabel('Degrees (0~360)'); grid on; xlim([0 t_end_days]);

%% 7. [부록] 순수 단주기 진동 확인 (Osculating - Mean)
figure('Name', 'Pure Short-Periodic Oscillations', 'Color', 'w');
sgtitle('Pure Short-Periodic Variations (Delta)', 'FontSize', 14, 'FontWeight', 'bold');

subplot(2, 2, 1);
plot(t_vec_days, a_osc_hist - a_mean_hist, 'b');
title('\Delta a (m)'); xlabel('Time (days)'); ylabel('m'); grid on;

subplot(2, 2, 2);
plot(t_vec_days, rad2deg(i_osc_hist - i_mean_hist), 'b');
title('\Delta i (deg)'); xlabel('Time (days)'); ylabel('deg'); grid on;

subplot(2, 2, 3);
plot(t_vec_days, rad2deg(RAAN_osc_hist - RAAN_mean_hist), 'b');
title('\Delta \Omega (deg)'); xlabel('Time (days)'); ylabel('deg'); grid on;

subplot(2, 2, 4);
plot(t_vec_days, rad2deg(u_osc_hist - u_mean_hist), 'b');
title('\Delta u (deg)'); xlabel('Time (days)'); ylabel('deg'); grid on;