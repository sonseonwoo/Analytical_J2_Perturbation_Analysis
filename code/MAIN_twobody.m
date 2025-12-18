%% Simple Two-Body Propagator

clc; clear all; close all; format longG;
global C
C = STEP1_constants();

%% 1. 초기 상수 정의
mu = 3.986004418e14; % 지구 중력 상수 (m^3/s^2)

%% 2. 초기 궤도 요소 (과제 조건)
a0_km = 6878;
e0 = 0.01;
i0_deg = 40;
RAAN0_deg = 120;
AoP0_deg = 20;
M0_deg = 60;

% --- 단위를 (m, rad)로 변환 ---
a0 = a0_km * 1000;
i0 = deg2rad(i0_deg);
RAAN0 = deg2rad(RAAN0_deg);
AoP0 = deg2rad(AoP0_deg);
M0 = deg2rad(M0_deg);


% --- M0 -> nu0 변환 ---
nu0 = M2nu(M0, e0,i0); % M2nu 함수가 필요합니다.

% --- COE -> 초기 상태 벡터 (r0, v0) 변환 ---
% coe2rv 함수가 필요합니다 (단위: m, m/s).
[r0_km, v0_km_s] = coe2rv(a0_km*(1-e0^2), e0, i0, RAAN0, AoP0, nu0, 0, 0, 0); % mu를 전달해야 할 수 있음
Y0 = [r0_km*1e3; v0_km_s*1e3];

%% 3. 전파 시간 설정 (예: 1일)
t_end_days = 1;
t_end_sec = t_end_days * 86400;
Step = 60; % 1분 간격
t_span = (0:Step:t_end_sec)'; % ode45는 t_span 전체를 입력으로 받음

%% 4. ODE 함수 정의 (파일 내부에 간단히 정의)
% 이 함수는 오직 2체 가속도만 계산합니다.
twoBodyODE = @(t, y) [ y(4:6); -mu / (norm(y(1:3))^3) * y(1:3) ];

%% 5. 수치 적분 수행 (ode45 사용)
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12); % 이전과 동일한 정밀도
fprintf('Starting simple Two-Body propagation...\n');
[t, yout] = ode89(twoBodyODE, t_span, Y0, options);
fprintf('Propagation finished.\n');

%% 6. 결과 분석: 에너지 계산 및 플롯
r_vectors_m = yout(:, 1:3);
v_vectors_m_s = yout(:, 4:6);
r_mag_m = sqrt(sum(r_vectors_m.^2, 2));
v_mag_m_s = sqrt(sum(v_vectors_m_s.^2, 2));
xi_m2_s2 = 0.5 * v_mag_m_s.^2 - mu ./ r_mag_m;
t_days = t / 86400;

figure('Name', 'Simple Two-Body: Specific Mechanical Energy');
plot(t_days, xi_m2_s2);
title('Simple Two-Body: Specific Mechanical Energy (\xi) vs. Time');
xlabel('Time (days)');
ylabel('Energy \xi (m^2/s^2)');
grid on;

% 참고: 에너지 변화를 더 잘 보기 위해 y축 범위 조정
mean_xi = mean(xi_m2_s2);
max_dev = max(abs(xi_m2_s2 - mean_xi));
if max_dev > 0 % 변화가 있을 경우에만 범위 조정
    ylim([mean_xi - 1.5*max_dev, mean_xi + 1.5*max_dev]);
end
fprintf('Mean Energy: %.10e\n', mean_xi);
fprintf('Max Energy Deviation: %.10e (%.4e %%)\n', max_dev, (max_dev/abs(mean_xi))*100);

%% 7. 결과 분석: 궤도 요소 계산 및 플롯 (선택 사항)
% rv2coe 함수가 필요합니다.
Eph_simple = zeros(length(t), 9); % 결과 저장 배열
for k = 1:length(t)
    [~, Eph_simple(k,1), Eph_simple(k,2), Eph_simple(k,3), Eph_simple(k,4), Eph_simple(k,5), Eph_simple(k,6), Eph_simple(k,7), ~, ~, ~] = ...
        rv2coe(yout(k,1:3)'/1000, yout(k,4:6)'/1000); % km, km/s 단위로 입력
end

figure('Name', 'Simple Two-Body: Orbital Elements', 'Color', 'w');
sgtitle('Simple Two-Body Propagation Results', 'FontSize', 14, 'FontWeight', 'bold');
subplot(3,2,1); plot(t_days, Eph_simple(:,1)); title('Semi-major Axis (a)'); xlabel('Time (days)'); ylabel('km'); grid on;
subplot(3,2,2); plot(t_days, Eph_simple(:,2)); title('Eccentricity (e)'); xlabel('Time (days)'); ylabel('Unitless'); grid on;
subplot(3,2,3); plot(t_days, rad2deg(Eph_simple(:,3))); title('Inclination (i)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;
subplot(3,2,4); plot(t_days, rad2deg(Eph_simple(:,4))); title('RAAN (\Omega)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;
subplot(3,2,5); plot(t_days, rad2deg(Eph_simple(:,5))); title('Argument of Perigee (\omega)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;
subplot(3,2,6); plot(t_days, rad2deg(Eph_simple(:,6))); title('true Anomaly (M)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;
subplot(3,2,6); plot(t_days, rad2deg(Eph_simple(:,7))); title('Mean Anomaly (M)'); xlabel('Time (days)'); ylabel('Degrees'); grid on;