%% 1. 초기 상수 정의
clc; clear all; close all; format longG;

% 지구 상수 (SI 단위, m/kg/s)
mu = 3.986004418e14; % (m^3/s^2)
R_e = 6378.137e3;    % 지구 반경 (m)
J2 = 0.00108263;     % J2 계수

%% 2. 초기 궤도 요소 (과제 조건 )
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

%% 3. 영년 변화율(Secular Rates) 계산
% 의 방정식을 사용합니다.

% --- 중간 변수 계산 ---
p0 = a0 * (1 - e0^2);         % Semi-latus rectum (m)
n0 = sqrt(mu / a0^3);         % 평균 운동 (rad/s) 의 n

% --- 궤도 요소별 영년 변화율 (rad/s) ---
a_rate = 0; % [cite: 242]
e_rate = 0; % [cite: 243]
i_rate = 0; % [cite: 244]

% RAAN 변화율 
RAAN_rate = -(3/2) * n0 * J2 * (R_e / p0)^2 * cos(i0);

% AoP 변화율 
AoP_rate = (3/4) * n0 * J2 * (R_e / p0)^2 * (5 * cos(i0)^2 - 1);

% 평균 근점 이각(M)의 변화율 
% M_dot = n_bar = n + <M_dot_0>
M_perturb_rate = (3/4) * n0 * J2 * (R_e / p0)^2 * sqrt(1 - e0^2) * (3 * cos(i0)^2 - 1);
M_rate_total = n0 + M_perturb_rate;

%% 4. 궤도 전파 (6개월)
% 6개월  (약 180일)간의 시간 벡터 생성
t_end_days = 1;
t_end_sec = t_end_days * 86400;
dt_sec = 60; % 1분 간격
t_vec = (0:dt_sec:t_end_sec)';

% --- 결과 저장을 위한 배열 초기화 ---
a_hist = ones(size(t_vec)) * a0;
e_hist = ones(size(t_vec)) * e0;
i_hist = ones(size(t_vec)) * i0;

% --- 전파 수행: Element = Element_0 + Rate * dt ---
RAAN_hist = RAAN0 + RAAN_rate * t_vec;
AoP_hist = AoP0 + AoP_rate * t_vec;
M_hist = M0 + M_rate_total * t_vec;
M0_hist = M0+M_perturb_rate*t_vec;

%% 5. 결과 시각화
t_vec_days = t_vec / 86400; % x축을 '일(days)'로 변경

figure('Name', 'Secular Perturbations (Analytical)');
sgtitle('J2 Secular Variations (Analytical Propagation)');

% --- RAAN ---
subplot(4, 1, 1);
plot(t_vec_days, rad2deg(RAAN_hist));
title('RAAN (\Omega) Evolution');
xlabel('Time (days)');
ylabel('Degrees');
grid on;

% --- Argument of Perigee ---
subplot(4, 1, 2);
plot(t_vec_days, rad2deg(AoP_hist));
title('Argument of Perigee (\omega) Evolution');
xlabel('Time (days)');
ylabel('Degrees');
grid on;

% --- Mean Anomaly (Wrapped) ---
subplot(4, 1, 3);
plot(t_vec_days, rad2deg(M_hist)); % 0~360도로 래핑
title('Mean Anomaly (M) Evolution (wrapped)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;

% --- Mean Anomaly (Wrapped) ---
subplot(4, 1, 4);
plot(t_vec_days, rad2deg(M0_hist)); % 0~360도로 래핑
title('Mean Anomaly (M0) Evolution (wrapped)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;



% --- 로그 출력: 6개월 후 변화량 ---
fprintf('--- 1day Secular Change (Analytical) ---\n');
fprintf('Initial RAAN: %.4f deg\n', RAAN0_deg);
fprintf('Final RAAN:   %.4f deg\n', rad2deg(RAAN_hist(end)));
fprintf('RAAN Rate:    %.4e deg/day\n\n', rad2deg(RAAN_rate) * 86400);

fprintf('Initial AoP: %.4f deg\n', AoP0_deg);
fprintf('Final AoP:   %.4f deg\n', rad2deg(AoP_hist(end)));
fprintf('AoP Rate:    %.4e deg/day\n\n', rad2deg(AoP_rate) * 86400);

% --- M0 출력 추가 ---
fprintf('Initial M0: %.4f deg\n', M0_deg);
fprintf('Final M0:   %.4f deg\n', rad2deg(M0_hist(end))); % 최종값은 unwrapped로 출력
fprintf('M0 Rate:    %.4e deg/day\n\n', rad2deg(M_perturb_rate) * 86400);

% --- M 출력 추가 ---
fprintf('Initial M: %.4f deg\n', M0_deg); % 초기 M은 M0와 같음
fprintf('Final M:   %.4f deg\n', rad2deg(M_hist(end))); % 최종값은 unwrapped로 출력
fprintf('M Rate (Total = n + M0_rate): %.4e deg/day\n\n', rad2deg(M_rate_total) * 86400);