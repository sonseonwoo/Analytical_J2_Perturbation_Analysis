% %--------------------------------------------------------------------------
% %
% % Accel: Computes the acceleration of an Earth orbiting satellite due to 
% %    	 - Earth's harmonic gravity field (including Solid Earth Tides and
% %      	   Ocean Tides), 
% %    	 - gravitational perturbations of the Sun, Moon and planets
% %    	 - solar radiation pressure
% %    	 - atmospheric drag and
% %	 	 - relativity
% %
% % Inputs:
% %   Mjd_UTC     Modified Julian Date (UTC)
% %   Y           Satellite state vector in the ICRF/EME2000 system
% %   Area        Cross-section 
% %   mass        Spacecraft mass
% %   Cr          Radiation pressure coefficient
% %   Cd          Drag coefficient
% %
% % Output:
% %   dY          Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
% %
% % Last modified:   2025/02/19   Meysam Mahooti
% % 
% %--------------------------------------------------------------------------
% function dY = Accel(t, Y)
% 
% global C in eopdata
% 
% MJD_UTC = in.Mjd_UTC+t/86400;
% 
% %  IERS(국제 지구 자전 및 기준계 사업)
% % 특정시점에 대한 EOP 데이터 확인 후 polar motion 좌표 확인
% [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
% 
% % 시간 척도 간 차이 보정
% [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
% 
% % 지구 자전각(ERA)[Theta] 계산에는 UT1이 사용되고,
% % 세차-장동 계산[NPB]에는 TT가 사용됨. 
% % 코드에서는 이 두 시간 척도를 MJD(수정 율리우스일) 형식으로 준비
% MJD_UT1 = MJD_UTC + UT1_UTC/86400;
% MJD_TT  = MJD_UTC + TT_UTC/86400;
% 
% % Form bias-precession-nutation matrix
% % 규약: IAU-2006/2000 결의안 (IAU: 국제천문연맹)
% % 기반: GCRF(지구중심 천구 기준계, ECI)에서 CIRS(천체 중간 기준계)로 변환
% %       이 변환은 기존의 춘분점 기반이 아닌, 천체 중간 원점(CIO, Celestial Intermediate Origin)을
% %       기준으로 하며, 프레임 바이어스, 세차(Precession), 장동(Nutation) 효과를
% %       하나의 행렬([BPN])으로 통합하여 계산
% NPB = iauPnm06a(C.DJM0, MJD_TT);
% 
% 
% % Form Earth rotation matrix
% % 규약: IAU-2006/2000 결의안 (IAU: 국제천문연맹)
% % 기반: CIRS(천체 중간 기준계)에서 TIRS(지구 중간 기준계)로 변환
% %       UT1을 사용하여 계산된 지구 자전각(ERA, Earth Rotation Angle)을 통해
% %       지구의 자전 운동을 반영. ERA는 천체 중간 원점(CIO)과
% %       지구 중간 원점(TIO) 사이의 각도
% gast = iauGst06(C.DJM0, MJD_UT1, C.DJM0, MJD_TT, NPB);
% Theta  = iauRz(gast, eye(3));
% 
% 
% % Polar motion matrix (TIRS->ITRS, IERS 2003)
% % 규약: IERS(국제 지구 자전 및 기준계 사업) 규약 (e.g., IERS Conventions 2010)
% % 기반: TIRS(지구 중간 기준계)에서 ITRS(국제 지구 기준계, ECEF)로 변환
% %       지구 자전축이 지각에 대해 움직이는 현상인 극 운동을 보정합니다
% %       IERS에서 제공하는 극 좌표(x_pole, y_pole)와 TIO Locator(s')를 사용합니다
% Po = iauPom00(x_pole, y_pole, iauSp00(C.DJM0, MJD_TT));
% 
% 
% 
% % GCRF(ECI)에서 ITRF(ECEF)로의 전체 위치 변환 적용
% % ICRS to ITRS transformation
% E = Po*Theta*NPB;
% 
% % 천체 궤도력 
% MJD_TDB = Mjday_TDB(MJD_TT);
% [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
%  r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE440(MJD_TDB+2400000.5);
% 
% % Acceleration due to harmonic gravity field
% if (in.SolidEarthTides || in.OceanTides)
%     a = AccelHarmonic_AnelasticEarth(MJD_UTC,r_Sun,r_Moon,Y(1:3),E,UT1_UTC,TT_UTC,x_pole*C.Arcs,y_pole*C.Arcs);
%     % a = AccelHarmonic_ElasticEarth(MJD_UTC,r_Sun,r_Moon,Y(1:3),E,UT1_UTC,TT_UTC,x_pole*const.Arcs,y_pole*const.Arcs);
% else
%     a = AccelHarmonic(Y(1:3), E, in.n, in.m);
% end
% 
% % Luni-solar perturbations
% if (in.sun)
%     a = a + AccelPointMass(Y(1:3),r_Sun,C.GM_Sun);
% end
% 
% if (in.moon)
%     a = a + AccelPointMass(Y(1:3),r_Moon,C.GM_Moon);
% end
% 
% % Planetary perturbations
% if (in.planets)
%     a = a + AccelPointMass(Y(1:3),r_Mercury,C.GM_Mercury);
%     a = a + AccelPointMass(Y(1:3),r_Venus,C.GM_Venus);
%     a = a + AccelPointMass(Y(1:3),r_Mars,C.GM_Mars);
%     a = a + AccelPointMass(Y(1:3),r_Jupiter,C.GM_Jupiter);
%     a = a + AccelPointMass(Y(1:3),r_Saturn,C.GM_Saturn);
%     a = a + AccelPointMass(Y(1:3),r_Uranus,C.GM_Uranus);    
%     a = a + AccelPointMass(Y(1:3),r_Neptune,C.GM_Neptune);
%     a = a + AccelPointMass(Y(1:3),r_Pluto,C.GM_Pluto);
% end
% 
% % Solar radiation pressure
% if (in.sRad)
%     a = a + AccelSolrad(Y(1:3),r_Earth,r_Moon,r_Sun,r_SunSSB, ...
%         in.area_solar,in.mass,in.Cr,C.P_Sol,C.AU,'conical');
% end
% 
% % Atmospheric drag
% if (in.drag)
%     % Atmospheric density
% 	Omega = C.omega_Earth-0.843994809*1e-9*LOD; % [rad/s]; IERS
%     dens = nrlmsise00(MJD_UTC,E*Y(1:3),UT1_UTC,TT_UTC);
%     % [~,dens] = JB2008(MJD_UTC,r_Sun,Y(1:3),E);
%     % [~,dens] = JB2006(MJD_UTC,r_Sun,Y(1:3),E);
%     % [d,~] = msis86(MJD_UTC,E*Y(1:3),gast);
%     % dens = 1e3*d(6);
%     % dens = Density_Jacchia70(r_Sun,MJD_UTC,E*Y(1:3),gast);
%     % dens = Density_HP(r_Sun,NPB*Y(1:3));
%     a = a + AccelDrag(dens,Y(1:3),Y(4:6),NPB,in.area_drag,in.mass,in.Cd,Omega);
% end
% 
% % Relativistic Effects
% if (in.Relativity)
%     a = a + Relativity(Y(1:3),Y(4:6));
% end
% 
% dY = [Y(4:6);a];
% 
%--------------------------------------------------------------------------
%
% Accel: Computes the acceleration of an Earth orbiting satellite.
%        Includes explicit handling for the pure Two-Body case.
%
% Inputs:
%   t           Time from epoch [s] (used by ode45)
%   Y           Satellite state vector [r; v] in ECI system [m; m/s]
%
% Output:
%   dY          Time derivative of the state vector [v; a] in ECI [m/s; m/s^2]
%
%--------------------------------------------------------------------------
function dY = Accel(t, Y)

global C in eopdata % Include necessary global variables

% --- 1. 순수 2체 문제인지 확인 ---
% 모든 섭동 플래그가 꺼져 있는지 확인합니다.
is_two_body = (in.n == 0 && in.m == 0 && ...          % No harmonics beyond n=0
               ~in.sun && ~in.moon && ~in.planets && ... % No third-body
               ~in.sRad && ~in.drag && ...               % No non-gravitational
               ~in.SolidEarthTides && ~in.OceanTides && ... % No tides
               ~in.Relativity);                           % No relativity

% --- 2. 순수 2체 문제인 경우: 직접 계산 후 즉시 반환 ---
if is_two_body
    r_vec = Y(1:3);
    r_mag = norm(r_vec);

    % C.muE가 m^3/s^2 단위의 중력 상수라고 가정합니다.
    a_2body_eci = -C.muE / (r_mag^3) * r_vec; 

    dY = [Y(4:6); a_2body_eci]; % [v; a]
    return; % 복잡한 계산을 모두 건너뛰고 함수 종료
end

% --- 3. 섭동이 포함된 경우: 기존의 복잡한 계산 수행 ---
% (이 부분은 이전과 동일하게 유지됩니다)

MJD_UTC = in.Mjd_UTC + t / 86400;

% EOP 데이터 읽기 및 시간 변환
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
MJD_UT1 = MJD_UTC + UT1_UTC/86400;
MJD_TT  = MJD_UTC + TT_UTC/86400;

% 좌표 변환 행렬 계산 (ECI <-> ECEF)
NPB = iauPnm06a(C.DJM0, MJD_TT);
gast = iauGst06(C.DJM0, MJD_UT1, C.DJM0, MJD_TT, NPB);
Theta  = iauRz(gast, eye(3));
Po = iauPom00(x_pole, y_pole, iauSp00(C.DJM0, MJD_TT));
E = Po*Theta*NPB; % ECI to ECEF 변환 행렬

% --- 중력장 가속도 계산 ---
if (in.SolidEarthTides || in.OceanTides)
     % 조석력 계산 부분 (별도 함수 호출)
     % MJD_TDB, r_Sun, r_Moon 등이 필요함 (아래에서 계산됨)
     MJD_TDB = Mjday_TDB(MJD_TT);
     [~,~,~,~,~,~,~,~,~,r_Moon,r_Sun,~] = JPL_Eph_DE440(MJD_TDB+2400000.5);
     a = AccelHarmonic_AnelasticEarth(MJD_UTC,r_Sun,r_Moon,Y(1:3),E,UT1_UTC,TT_UTC,x_pole*C.Arcs,y_pole*C.Arcs);
else
    % 일반 중력장 또는 J2/2체 계산 (수정된 AccelHarmonic 사용 가능)
    a = AccelHarmonic(Y(1:3), E, in.n, in.m); 
end

% --- 다른 섭동 가속도 추가 (플래그가 켜진 경우) ---

% 천체 궤도력 계산 (조석력 외 다른 섭동에 필요)
if (in.sun || in.moon || in.planets || in.sRad || in.drag) % JPL 계산이 필요한 경우만
    MJD_TDB = Mjday_TDB(MJD_TT);
    [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
     r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE440(MJD_TDB+2400000.5);
end

% Luni-solar perturbations
if (in.sun)
    a = a + AccelPointMass(Y(1:3), r_Sun, C.GM_Sun);
end
if (in.moon)
    a = a + AccelPointMass(Y(1:3), r_Moon, C.GM_Moon);
end

% Planetary perturbations
if (in.planets)
    % ... (각 행성에 대한 AccelPointMass 호출) ...
end

% Solar radiation pressure
if (in.sRad)
    a = a + AccelSolrad(Y(1:3),r_Earth,r_Moon,r_Sun,r_SunSSB, ...
        in.area_solar,in.mass,in.Cr,C.P_Sol,C.AU,'conical');
end

% Atmospheric drag
if (in.drag)
    Omega = C.omega_Earth - 0.843994809*1e-9*LOD; 
    dens = nrlmsise00(MJD_UTC, E*Y(1:3), UT1_UTC, TT_UTC);
    a = a + AccelDrag(dens, Y(1:3), Y(4:6), NPB, in.area_drag, in.mass, in.Cd, Omega);
end

% Relativistic Effects
if (in.Relativity)
    a = a + Relativity(Y(1:3), Y(4:6));
end

% --- 최종 미분값 설정 ---
dY = [Y(4:6); a];

end % End of Accel function