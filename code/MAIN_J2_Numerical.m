
%%
clc
clear all
format longG
close all
global C in Cnm Snm eopdata swdata SOLdata DTCdata APdata PC
global_variable;
C = STEP1_constants();
% ─── 0. propagator environment & dynamic modeling setting ─────────
in.T = datetime('2025-04-28 00:00:00','TimeZone','UTC'); % 관측 요청 시각
in.TV = datevec(in.T);
in.Mjd_UTC =Mjday(in.TV(1), in.TV(2),in.TV(3), in.TV(4),in.TV(5),in.TV(6)); 
in.n       = 2; % J2 섭동만
in.m       = 0;
in.sun     = 0; in.moon = 0; in.planets = 0;
in.sRad    = 0; in.drag  = 0;
in.SolidEarthTides = 0; in.OceanTides = 0; in.Relativity = 0;
in.mass = 10; 
in.options = odeset('RelTol',1e-5,'AbsTol',1e-5);

%% 1. 초기 궤도 요소 설정
a = 6878; e = 0.01;
inc = deg2rad(40); RAAN = deg2rad(120);
AoP =deg2rad(20); Mo = deg2rad(60);
nu = M2nu(Mo,e,inc); % inc 인자 제거
p = a*(1-e^2);

%% 2. 수치적 궤도 전파 (Cowell) - 30일
[r0_ijk,v0_ijk]  = coe2rv(p,e,inc,RAAN,AoP,nu,0,0,0); 
Y0 = [r0_ijk;v0_ijk].*1e3; % km -> m

% 3.2절 요구사항: 30일
Step   = 100;  
N_Step = 1*864; % (30 * 86400 / 100)    
t_vec = (0:Step:N_Step*Step)'; 

fprintf('Starting ODE45 propagation for 000 days...\n');
tic
[t,yout] = ode89(@Accel, t_vec, Y0, in.options); 
fprintf('Propagation finished.\n');
toc
tic
% 상태벡터 -> 순간 접촉 궤도 요소 (Osculating) 변환
[~,~, Eph_osc] = processOrbitData(yout, t);
% Eph_osc 열: [~, ~, a(m), ecc, incl(rad), RAAN(rad), AoP(rad), nu(rad), M(rad)]
t_days = t / 86400; % x축 (일)
toc
%% 3.1 영년 변화 (Secular Trend) 분석
% --- 모든 궤도요소의 영년 추세선(평균값) 계산 ---

% a (장반경)
a_osc_m = Eph_osc(:,3);
p_a = polyfit(t, a_osc_m, 1);
a_trend_m = polyval(p_a, t);

% e (이심률)
e_osc = Eph_osc(:,4);
p_e = polyfit(t, e_osc, 1);
e_trend = polyval(p_e, t);

% i (궤도 경사각)
i_osc_rad = Eph_osc(:,5);
p_i = polyfit(t, i_osc_rad, 1);
i_trend_rad = polyval(p_i, t);

% RAAN
RAAN_osc_rad = Eph_osc(:,6);
RAAN_unwrapped_rad = unwrap(RAAN_osc_rad);
p_RAAN = polyfit(t, RAAN_unwrapped_rad, 1);
RAAN_trend_rad = polyval(p_RAAN, t);

% AoP
AoP_osc_rad = Eph_osc(:,7);
AoP_unwrapped_rad = unwrap(AoP_osc_rad);
p_AoP = polyfit(t, AoP_unwrapped_rad, 1);
AoP_trend_rad = polyval(p_AoP, t);

% nu (진근점 이각) --- [추가됨]
nu_osc_rad = Eph_osc(:,8);
nu_unwrapped_rad = unwrap(nu_osc_rad);
p_nu = polyfit(t, nu_unwrapped_rad, 1);
nu_trend_rad = polyval(p_nu, t);

% M (평균 근점 이각) --- [추가됨]
M_osc_rad = Eph_osc(:,9);
M_unwrapped_rad = unwrap(M_osc_rad);
p_M = polyfit(t, M_unwrapped_rad, 1);
M_trend_rad = polyval(p_M, t);


% M0 (평균 근점 이각 초기값) 
M_trend2 = polyval([sqrt(C.muE*1e-9/p_a(2)^3) 0],t);
M0_unwrapped_rad = M_unwrapped_rad-M_trend2;
p_M0 = polyfit(t,M0_unwrapped_rad,1);
M0_trend_rad = polyval(p_M0,t);

% Xi
Xi_osc_j2 = Eph_osc(:,14);

Xi_osc = Eph_osc(:,13);



%%

% 2. 궤도 요소 플로팅
figure('Name', 'Osculating Orbital Elements (J2 Perturbation)', 'Color', 'w');
sgtitle('Osculating Orbital Elements (J2 Perturbation Only)', 'FontSize', 14, 'FontWeight', 'bold');

% --- Subplot 1: Semi-major Axis (a) ---
subplot(4,2,1)
plot(t_days, Eph_osc(:,3) ); % m -> km
title('Semi-major Axis (a)');
xlabel('Time (days)');
ylabel('km');
grid on;

% --- Subplot 2: Eccentricity (e) ---
subplot(4,2,2)
plot(t_days, Eph_osc(:,4));
title('Eccentricity (e)');
xlabel('Time (days)');
ylabel('Unitless');
grid on;

% --- Subplot 3: Inclination (i) ---
subplot(4,2,3)
plot(t_days, rad2deg(Eph_osc(:,5))); % rad -> deg
title('Inclination (i)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;

% --- Subplot 4: RAAN (Omega) ---
subplot(4,2,4)
plot(t_days, rad2deg(RAAN_unwrapped_rad)); % rad -> deg
title('RAAN (\Omega)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;

% --- Subplot 5: Argument of Perigee (AoP) ---
subplot(4,2,5)
plot(t_days, rad2deg(AoP_unwrapped_rad)); % rad -> deg
title('Argument of Perigee (\omega)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;

% --- Subplot 6: True Anomaly (nu) ---
subplot(4,2,6)
plot(t_days, rad2deg(nu_unwrapped_rad)); % rad -> deg
title('True Anomaly (\nu)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;

% --- Subplot 7: Mean Anomaly (M) ---
subplot(4,2,7)
plot(t_days, rad2deg(M_unwrapped_rad)); % rad -> deg
title('Mean Anomaly (M)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;

% --- Subplot 8: Mean Anomaly (M0) ---
subplot(4,2,8)
plot(t_days, rad2deg(M0_unwrapped_rad)); % rad -> deg
title('Initial Mean Anomaly (M0)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;




%%
% --- 3.1. 피팅 결과 시각화 (진동 vs 추세) ---
figure('Name', 'Numerical Trend vs. Osculating Data (3.1)');
sgtitle('Osculating vs. Secular Trend (90 Days)', 'FontSize', 14);

subplot(4,2,1)
plot(t_days, a_osc_m , 'Color', [0.7 0.7 0.7]); hold on;
plot(t_days, a_trend_m , 'r-', 'LineWidth', 1.5);
title('Semi-major Axis (a)'); xlabel('Time (days)'); ylabel('km');
legend('Osculating', 'Secular Trend', 'Location', 'best'); grid on;

subplot(4,2,2)
plot(t_days, e_osc, 'Color', [0.7 0.7 0.7]); hold on;
plot(t_days, e_trend, 'r-', 'LineWidth', 1.5);
title('Eccentricity (e)'); xlabel('Time (days)'); ylabel('Unitless');
grid on;

subplot(4,2,3)
plot(t_days, rad2deg(i_osc_rad), 'Color', [0.7 0.7 0.7]); hold on;
plot(t_days, rad2deg(i_trend_rad), 'r-', 'LineWidth', 1.5);
title('Inclination (i)'); xlabel('Time (days)'); ylabel('Degrees');
grid on;

subplot(4,2,4)
plot(t_days, rad2deg(RAAN_unwrapped_rad), 'Color', [0.7 0.7 0.7]); hold on;
plot(t_days, rad2deg(RAAN_trend_rad), 'r-', 'LineWidth', 1.5);
title('RAAN (\Omega)'); xlabel('Time (days)'); ylabel('Degrees');
grid on;

subplot(4,2,5)
plot(t_days, rad2deg(AoP_unwrapped_rad), 'Color', [0.7 0.7 0.7]); hold on;
plot(t_days, rad2deg(AoP_trend_rad), 'r-', 'LineWidth', 1.5);
title('Argument of Perigee (\omega)'); xlabel('Time (days)'); ylabel('Degrees');
grid on;

% --- M (평균 근점 이각) 플롯 [추가됨] ---
subplot(4,2,6)
plot(t_days, rad2deg(M_unwrapped_rad), 'Color', [0.7 0.7 0.7]); hold on;
plot(t_days, rad2deg(M_trend_rad), 'r-', 'LineWidth', 1.5);
title('Mean Anomaly (M)'); xlabel('Time (days)'); ylabel('Degrees');
grid on;

% --- M0 (평균 근점 이각) 플롯 [추가됨] ---
subplot(4,2,7)
plot(t_days, rad2deg(M0_unwrapped_rad), 'Color', [0.7 0.7 0.7]); hold on;
plot(t_days, rad2deg(M0_trend_rad), 'r-', 'LineWidth', 1.5);
title('Mean Anomaly (M0)'); xlabel('Time (days)'); ylabel('Degrees');
grid on;




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.2 장주기 변화 (수치적 분석) - "추세 제거" 플롯
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- "순간값(Osculating) - 추세선(Trend)" 계산 ---
a_detrended_m = a_osc_m - a_trend_m;
e_detrended = e_osc - e_trend;
i_detrended_rad = i_osc_rad - i_trend_rad;
RAAN_detrended_rad = RAAN_unwrapped_rad - RAAN_trend_rad;
AoP_detrended_rad = AoP_unwrapped_rad - AoP_trend_rad;
M_detrended_rad = M_unwrapped_rad - M_trend_rad; % [추가됨]
M0_detrended_rad = M0_unwrapped_rad -M0_trend_rad;


% --- 3.2. "추세 제거" 플롯 ---
figure('Name', 'Detrended Elements (3.2 Long-Periodic Analysis)', 'Color', 'w');
sgtitle('Osculating Elements minus Secular Trend (Pure Periodic Variations)', 'FontSize', 14, 'FontWeight', 'bold');

% --- Subplot 1: Semi-major Axis (a) ---
subplot(4,2,1)
plot(t_days, a_detrended_m ); % m -> km
title('Detrended Semi-major Axis (a)');
xlabel('Time (days)');
ylabel('km');
grid on;

% --- Subplot 2: Eccentricity (e) ---
subplot(4,2,2)
plot(t_days, e_detrended);
title('Detrended Eccentricity (e)');
xlabel('Time (days)');
ylabel('Unitless');
grid on;

% --- Subplot 3: Inclination (i) ---
subplot(4,2,3)
plot(t_days, rad2deg(i_detrended_rad)); % rad -> deg
title('Detrended Inclination (i)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;

% --- Subplot 4: RAAN (Omega) ---
subplot(4,2,4)
plot(t_days, rad2deg(RAAN_detrended_rad)); % rad -> deg
title('Detrended RAAN (\Omega)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;

% --- Subplot 5: Argument of Perigee (AoP) ---
subplot(4,2,5)
plot(t_days, rad2deg(AoP_detrended_rad)); % rad -> deg
title('Detrended Argument of Perigee (\omega)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;

% --- Subplot 6: Mean Anomaly (M) [추가됨] ---
subplot(4,2,6)
plot(t_days, rad2deg(M_detrended_rad)); % rad -> deg

title('Detrended Mean Anomaly (M)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;

% --- Subplot 6: Mean Anomaly (M0) [추가됨] ---
subplot(4,2,7)
plot(t_days, rad2deg(M0_detrended_rad)); % rad -> deg

title('Detrended Mean Anomaly (M0)');
xlabel('Time (days)');
ylabel('Degrees');
grid on;




%% 4. 수치적 결과 출력 (30일 기준)
% (기존 3.1절에서 계산한 값들을 사용)

% --- 속도(rate) 계산 ---
RAAN_rate_num_deg_day = rad2deg(p_RAAN(1)) * 86400;
AoP_rate_num_deg_day = rad2deg(p_AoP(1)) * 86400;
M_rate_num_deg_day = rad2deg(p_M(1)) * 86400; % [추가됨]
M0_rate_num_deg_day = rad2deg(p_M0(1))*86400;

% --- 초기/최종 값 (추세선 기준) ---
RAAN0_num_deg = rad2deg(RAAN_trend_rad(1));
RAAN_final_num_deg = rad2deg(RAAN_trend_rad(end));
AoP0_num_deg = rad2deg(AoP_trend_rad(1));
AoP_final_num_deg = rad2deg(AoP_trend_rad(end));
M_num_deg = rad2deg(M_trend_rad(1)); % [추가됨]
M_final_num_deg = rad2deg(M_trend_rad(end)); % [추가됨]
M0_num_deg = rad2deg(M0_trend_rad(1));
M0_final_num_deg = rad2deg(M0_trend_rad(end));


% --- 로그 출력 ---
fprintf('\n--- 30-Day Secular Change (Numerical/Cowell) ---\n');
fprintf('Initial RAAN (mean): %.4f deg\n', RAAN0_num_deg);
fprintf('Final RAAN (mean):   %.4f deg\n', RAAN_final_num_deg);
fprintf('RAAN Rate (mean):    %.4e deg/day\n\n', RAAN_rate_num_deg_day);

fprintf('Initial AoP (mean): %.4f deg\n', AoP0_num_deg);
fprintf('Final AoP (mean):   %.4f deg\n', AoP_final_num_deg);
fprintf('AoP Rate (mean):    %.4e deg/day\n\n', AoP_rate_num_deg_day);

fprintf('Initial M0 (mean): %.4f deg\n', M0_num_deg); % [추가됨]
fprintf('Final M0 (mean):   %.4f deg\n', M0_final_num_deg); % [추가됨]
fprintf('M0 Rate (mean):    %.4e deg/day\n\n', M0_rate_num_deg_day); % [추가됨]

fprintf('Initial M (mean): %.4f deg\n', M_num_deg); % [추가됨]
fprintf('Final M (mean):   %.4f deg\n', M_final_num_deg); % [추가됨]
fprintf('M Rate (mean):    %.4e deg/day\n\n', M_rate_num_deg_day); % [추가됨]
