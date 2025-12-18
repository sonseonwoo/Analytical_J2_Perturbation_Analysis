function [Eph_eci, Eph_ecef, Eph_osc] = processOrbitData(yout, t)
    % 궤도 전파 결과(yout, t)로부터 ECI, ECEF, Osculating Elements 데이터를 생성합니다.
    % Inputs:
    %   yout    : [N x 6] ECI 상태 벡터 (m, m/s)
    %   t       : [N x 1] 시간 벡터 (s)
    %   C       : 상수 구조체 (muE, RE, J2 등, m-s 단위)
    %   Mjd_UTC : 기준 시각 (Modified Julian Date)
    % Outputs:
    %   Eph_eci  : [N x 8] ECI 데이터 [t, r, v, a]
    %   Eph_ecef : [N x 8] ECEF 데이터 [t, r, v, alt]
    %   Eph_osc  : [N x 14] 궤도 요소 데이터 [t, p, a, ecc, ...]
    global C in
    Mjd_UTC = in.Mjd_UTC;%+((600+5837)/86400);
    [n_steps, ~] = size(yout);

    % --- 1. Eph_eci 계산 (장반경 포함) ---
    Eph_eci = zeros(n_steps, 8);
    for i = 1:n_steps
        r_vec = yout(i, 1:3);
        v_vec = yout(i, 4:6);
        r_mag = norm(r_vec);
        v_mag_sq = dot(v_vec, v_vec);
        a = 1 / ( (2/r_mag) - (v_mag_sq / C.muE) );
        Eph_eci(i, :) = [t(i), yout(i,:), a];
    end

    % --- 2. Eph_ecef 계산 (고도 포함) ---
    Eph_ecef = zeros(n_steps, 8);
    for i = 1:n_steps
        Eph_ecef(i,1) = Eph_eci(i,1);
        state_eci_m = Eph_eci(i, 2:7);
        state_ecef_m = ECI2ECEF(Mjd_UTC + Eph_ecef(i,1)/86400, state_eci_m);
        lla = ecef2lla(state_ecef_m(1:3)'); % ecef2lla는 보통 m 단위를 입력으로 받음
        Eph_ecef(i,2:7) = state_ecef_m;
        Eph_ecef(i,8) = lla(3); % 고도(m)
    end

    % --- 3. Eph_osc 계산 (궤도 요소) ---
    Eph_osc = zeros(n_steps, 14);
    for i = 1:n_steps
        r_vec_km = Eph_eci(i, 2:4) / 1000;
        v_vec_km = Eph_eci(i, 5:7) / 1000;
        
        [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper,xi,xi_J2] = ...
            rv2coe(r_vec_km, v_vec_km);
        
        Eph_osc(i, :) = [Eph_eci(i, 1), p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper, xi, xi_J2];
    end
end