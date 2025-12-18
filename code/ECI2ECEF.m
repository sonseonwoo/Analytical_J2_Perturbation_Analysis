%--------------------------------------------------------------------------
%
% ECI2ECEF: Transforms Earth Centered Inertial (ECI) coordinates to Earth
%           Centered Earth Fixed (ECEF) coordinates
%
% input : 
%         Y0   : %[1 x 6]
% Last modified:   2022/11/07   Meysam Mahooti
% 
%--------------------------------------------------------------------------
function Y = ECI2ECEF(MJD_UTC, Y0)



global C eopdata

%  IERS(국제 지구 자전 및 기준계 사업)
% 특정시점에 대한 EOP 데이터 확인 후 polar motion 좌표 확인
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');

% 시간 척도 간 차이 보정
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);

% 지구 자전각(ERA)[Theta] 계산에는 UT1이 사용되고,
% 세차-장동 계산[NPB]에는 TT가 사용됨. 
% 코드에서는 이 두 시간 척도를 MJD(수정 율리우스일) 형식으로 준비

MJD_UT1 = MJD_UTC + UT1_UTC/86400;
MJD_TT = MJD_UTC + TT_UTC/86400;

% Form bias-precession-nutation matrix
% 규약: IAU-2006/2000 결의안 (IAU: 국제천문연맹)
% 기반: GCRF(지구중심 천구 기준계, ECI)에서 CIRS(천체 중간 기준계)로 변환
%       이 변환은 기존의 춘분점 기반이 아닌, 천체 중간 원점(CIO, Celestial Intermediate Origin)을
%       기준으로 하며, 프레임 바이어스, 세차(Precession), 장동(Nutation) 효과를
%       하나의 행렬([BPN])으로 통합하여 계산
NPB = iauPnm06a(C.DJM0, MJD_TT);


% Form Earth rotation matrix
% 규약: IAU-2006/2000 결의안 (IAU: 국제천문연맹)
% 기반: CIRS(천체 중간 기준계)에서 TIRS(지구 중간 기준계)로 변환
%       UT1을 사용하여 계산된 지구 자전각(ERA, Earth Rotation Angle)을 통해
%       지구의 자전 운동을 반영. ERA는 천체 중간 원점(CIO)과
%       지구 중간 원점(TIO) 사이의 각도
Theta = iauRz( iauGst06(C.DJM0, MJD_UT1, C.DJM0, MJD_TT, NPB),eye(3) );


% Polar motion matrix (TIRS->ITRS, IERS 2003)
% 규약: IERS(국제 지구 자전 및 기준계 사업) 규약 (e.g., IERS Conventions 2010)
% 기반: TIRS(지구 중간 기준계)에서 ITRS(국제 지구 기준계, ECEF)로 변환
%       지구 자전축이 지각에 대해 움직이는 현상인 극 운동을 보정합니다
%       IERS에서 제공하는 극 좌표(x_pole, y_pole)와 TIO Locator(s')를 사용합니다
Po = iauPom00(x_pole, y_pole, iauSp00(C.DJM0, MJD_TT));




% ICRS to ITRS transformation matrix and derivative
S = zeros(3);
S(1,2) = 1;
S(2,1) = -1;
Omega  = C.omega_Earth-0.843994809*1e-9*LOD;     % [rad/s]; IERS
dTheta = Omega*S*Theta;           				 % Derivative of Earth rotation matrix [1/s]


% GCRF(ECI)에서 ITRF(ECEF)로의 전체 위치 변환 적용
U      = Po*Theta*NPB;                			 % ICRS to ITRS transformation
% GCRF(ECI)에서 ITRF(ECEF)로의 전체 속도 변환 적용
dU     = Po*dTheta*NPB;               			 % Derivative [1/s]

% Transformation from ICRS to WGS
r = U*Y0(1:3)';                                  % [3 x 3] [3 x 1]
v = U*Y0(4:6)' + dU*Y0(1:3)';                    % [3 x 3] [3 x 1]
Y = [r;v];                                       % [6 x 1]

