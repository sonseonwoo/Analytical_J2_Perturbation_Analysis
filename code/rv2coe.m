function [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper,sme,sme_J2] = rv2coe(r, v)
% r, v, mu를 입력으로 받아 궤도 요소를 계산하는 강건성이 보완된 함수
% Outputs are in SI units (km, rad, rad/s)
%
% author        : david vallado (original logic), user (refactored)
%
% inputs          description                    range / units
%   r           - ijk position vector            km
%   v           - ijk velocity vector            km / s
%   mu          - gravitational parameter        km^3 / s^2
%
% outputs       :
%   p           - semilatus rectum               km
%   a           - semimajor axis                 km
%   ecc         - eccentricity
%   incl        - inclination                    0.0 to pi rad
%   omega       - longitude of ascending node    0.0 to 2pi rad
%   argp        - argument of perigee            0.0 to 2pi rad
%   nu          - true anomaly                   0.0 to 2pi rad
%   m           - mean anomaly                   0.0 to 2pi rad
%   arglat      - argument of latitude      (ci) 0.0 to 2pi rad
%   truelon     - true longitude            (ce) 0.0 to 2pi rad
%   lonper      - longitude of periapsis    (ee) 0.0 to 2pi rad
%   xi          - specific mechanical energy     km^2 / s^2 (user requested)
%   xi_J2       - specific mech. energy w/ J2    km^2 / s^2 (user requested)
global C 
% --- 1. 초기 상수 및 변수 설정 ---
small     = 1e-3;      % 부동소수점 비교를 위한 작은 값
undefined = NaN;        % 정의되지 않은 값은 NaN으로 처리
twopi     = 2.0 * pi;


% J2 섭동과 관련된 상수
mu = C.muE/1e9;    %[km^3/s^2]
J2 = C.J2;
R_earth = C.RE/1e3;


% 출력 변수 초기화
p = undefined; a = undefined; ecc = undefined; incl = undefined;
omega = undefined; argp = undefined; nu = undefined; m = undefined;
arglat = undefined; truelon = undefined; lonper = undefined;


% --- 2. 기본 벡터 계산 ---
magr = norm(r);
magv = norm(v);
rdotv = dot(r, v);

hbar = cross(r, v); % 각운동량 벡터
magh = norm(hbar);

% 궤도가 존재할 수 있는 기본 조건 확인
if magh > small
    nbar = [-hbar(2), hbar(1), 0.0]; % 승교점 벡터
    magn = norm(nbar);

    ebar = ((magv^2 - mu/magr)*r - rdotv*v) / mu; % 이심률 벡터
    ecc = norm(ebar);

    % --- 3. 에너지 및 주요 궤도 요소 계산 ---
    sme = 0.5*magv^2 - mu/magr; % Specific mechanical energy (= sme)
    sme_J2 = sme + (J2*mu*R_earth^2)/(2*magr^3)*(1-3*(r(3)^2/magr^2));

    if abs(1.0 - ecc) > small % 포물선 궤도가 아닐 때
        a = -mu / (2.0 * sme);
        p = a * (1.0 - ecc^2);
    else % 포물선 궤도일 때
        a = inf;
        p = magh^2 / mu;
    end

    % --- 4. 궤도 경사각 (inclination) 계산 ---
    % acos의 입력값이 -1~1 범위를 벗어나지 않도록 clamp
    temp = hbar(3)/magh;
    if(abs(temp) > 1.0)
        temp = sign(temp);
    end
    incl = acos(temp);

    % --- 5. 특수 궤도 유형 분류 ---
    if ecc < small % 원궤도 계열
        if (incl < small) || (abs(incl - pi) < small)
            typeorbit = 'ce'; % Circular Equatorial
        else
            typeorbit = 'ci'; % Circular Inclined
        end
    else % 비원궤도 계열
        if (incl < small) || (abs(incl - pi) < small)
            typeorbit = 'ee'; % Elliptical Equatorial
        else
            typeorbit = 'ei'; % Elliptical Inclined
        end
    end

    % --- 6. 궤도 요소 계산 (궤도 유형에 따라) ---
    % 승교점 경도 (RAAN)
    if magn > small
        temp = nbar(1) / magn;
        if(abs(temp) > 1.0)
            temp = sign(temp);
        end
        omega = acos(temp);
        if nbar(2) < 0.0
            omega = twopi - omega;
        end
    end

    % 근지점 인수 (Argument of Perigee)
    if strcmp(typeorbit, 'ei')
        temp = dot(nbar, ebar) / (magn * ecc);
        if(abs(temp) > 1.0)
            temp = sign(temp);
        end
        argp = acos(temp);
        if ebar(3) < 0.0
            argp = twopi - argp;
        end
    end

    % 진근점 이각 (True Anomaly)
    if typeorbit(1) == 'e' % ee, ei
        temp = dot(ebar, r) / (ecc * magr);
        if(abs(temp) > 1.0)
            temp = sign(temp);
        end
        nu = acos(temp);
        if rdotv < 0.0
            nu = twopi - nu;
        end
    end
    
    % --- 7. 특수 궤도 각도 계산 ---
    % 궤도 위도 인수 (Argument of Latitude)
    if strcmp(typeorbit, 'ci')
        temp = dot(nbar, r) / (magn * magr);
        if(abs(temp) > 1.0)
            temp = sign(temp);
        end
        arglat = acos(temp);
        if r(3) < 0.0
            arglat = twopi - arglat;
        end
    end
    
    % 근지점 경도 (Longitude of Periapsis)
    if strcmp(typeorbit, 'ee') && ecc > small
        temp = ebar(1) / ecc;
        if(abs(temp) > 1.0)
            temp = sign(temp);
        end
        lonper = acos(temp);
        if ebar(2) < 0.0
            lonper = twopi - lonper;
        end
    end

    % 진경도 (True Longitude)
    if strcmp(typeorbit, 'ce') && magr > small
        temp = r(1) / magr;
        if(abs(temp) > 1.0)
            temp = sign(temp);
        end
        truelon = acos(temp);
        if r(2) < 0.0
            truelon = twopi - truelon;
        end
    end

    % --- 8. 평균 근점 이각 (Mean Anomaly) 계산 ---
    if typeorbit(1) == 'e' % ee, ei
        % True Anomaly -> Eccentric Anomaly -> Mean Anomaly
        E = 2 * atan(sqrt((1-ecc)/(1+ecc)) * tan(nu/2));
        m = E - ecc * sin(E);
        
        % m을 [0, 2*pi] 범위로 변환 (음수 각도를 양수 각도로 래핑)
        m = mod(m, 2*pi);
        
    else % ci, ce
        % 원궤도에서는 Mean Anomaly가 정의되지 않지만, 
        % 종종 Argument of Latitude나 True Longitude를 대신 사용
        if strcmp(typeorbit, 'ci')
            m = arglat;
        elseif strcmp(typeorbit, 'ce')
            m = truelon;
        end
    end
end
end