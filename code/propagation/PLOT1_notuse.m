tic
a          = best.osc.a*1e-3;
inc        = best.osc.i;
RAAN       = best.osc.RAAN;

if isfield(best.osc, 'u')
    % best.sol.u0가 존재할 때 실행
    arglat     = best.osc.u;
    % r0 [3 x 1] km v0 [3 x 1] km/s
    [r0_ijk,v0_ijk]  = coe2rv(a,0,inc,RAAN,0,0,arglat,0,0); 
else
    % 존재하지 않을 때 실행
    e = best.osc.e;
    AoP = best.osc.AoP0;
    nu = best.osc.nu0;
    [r0_ijk,v0_ijk]  = coe2rv(a,e,inc,RAAN,AoP,nu,0,0,0); 

end

% a          = best.sol.a*1e-3;
% inc        = best.sol.i;
% RAAN       = best.sol.RAAN0;
% 
% if isfield(best.sol, 'u0')
%     % best.sol.u0가 존재할 때 실행
%     arglat     = best.sol.u0;
%     % r0 [3 x 1] km v0 [3 x 1] km/s
%     [r0_ijk,v0_ijk]  = coe2rv(a,0,inc,RAAN,0,0,arglat,0,0); 
% else
%     % 존재하지 않을 때 실행
%     e = best.sol.e;
%     AoP = best.sol.AoP0;
%     nu = best.sol.nu0;
%     [r0_ijk,v0_ijk]  = coe2rv(a,e,inc,RAAN,AoP,nu,0,0,0); 
% 
% end

%%
% ─── 0. propagator environment & dynamic modeling setting ─────────
in.n       = 2;
in.m       = 0;
in.sun     = 0; in.moon = 0; in.planets = 0;
in.sRad    = 0; in.drag  = 0;
in.SolidEarthTides = 0; in.OceanTides = 0; in.Relativity = 0;


% Y0 [6 x 1] [m]
Y0 = [r0_ijk;v0_ijk].*1e3;
Step   = 1;   % [s] integration step size
% N_Step = 15*(2*pi/sqrt(C.muE/best.osc.a^3))+best.tp; % number of integration steps 
N_Step = best.tp-4; % number of integration steps 
[t,yout] = ode113(@Accel,(0:Step:N_Step*Step),Y0,in.options); 

[Eph_eci, Eph_ecef, Eph_osc] = processOrbitData(yout, t);

r_ecef_km = (Eph_ecef(:,2:4))./1e3;      % [N x 3]

% 
% GroundHist = ecef2lla(Eph_ecef(:,2:4));

% WGS84와 동일한 평균 반지름을 가진 구체 생성
earthSphere = referenceSphere('Earth','meter');
% Eph_ecef 행렬에서 x, y, z 좌표를 각각의 열 벡터로 분리
x = Eph_ecef(:, 2);
y = Eph_ecef(:, 3);
z = Eph_ecef(:, 4);
% x, y, z를 개별 인수로 전달하여 함수 호출
[lat, lon, alt] = ecef2geodetic(x, y, z, earthSphere);
% GroundHist 변수에 결과 저장 (필요한 경우)
GroundHist = [rad2deg(lat), rad2deg(lon), alt]; 


GroundHist(:,2) = mod(GroundHist(:,2),360);
toc
%%


X0lla = GroundHist(1,:);
Xflla = GroundHist(end,:);

tlat = rad2deg(in.target.lat);
tlon = rad2deg(in.target.lon);
llat = rad2deg(in.launch.lat);
llon = rad2deg(in.launch.lon);
[llatc, llonc] = scircle1(llat, llon, km2deg(2*C.RE*1e-3*abs(best.resMax(1))*cos(in.launch.lat_gc)));
[slatcC2, sloncC2] = scircle1(X0lla(1), X0lla(2), km2deg(2*C.RE*1e-3*abs(best.resMax(2))*cosd(X0lla(1))));
% [tlatcC1, tloncC1] = scircle1(tlat, tlon, km2deg(2*C.RE*1e-3*abs(best.resMax(1))*cos(in.target.lat_gc)));%1e9*
[tlatcC1, tloncC1] = scircle1(tlat, tlon, km2deg(2*C.RE*1e-3*abs(0.02)*cos(in.target.lat_gc)));%1e9*

[tlatcC3, tloncC3] = scircle1(tlat, tlon, km2deg(2*C.RE*1e-3*abs(best.resMax(2))*cos(in.target.lat_gc)));%1e9*


target_lladeg = [tlat, tlon, 459.22];
launch_lladeg = [llat, llon, 0];

[eta_deg, vis_horizon, vis_mission] = fov_nadir_conical(r_ecef_km, target_lladeg, C.RE/1e3, 1, true);
% [eta_deg2, vis_horizon2, vis_mission2] = fov_nadir_conical(r_ecef_km, launch_lladeg, C.RE/1e3, 10, true);
orb_elements = [a/1000, inc]; % a는 km, inc는 rad
target_lla_gt = [tlat, tlon];
launch_lla_gt = [rad2deg(in.launch.lat), rad2deg(in.launch.lon)];

target_lla_3d = [tlat, tlon, target_lladeg(3)]; % 고도 단위를 m로
launch_lla_3d = [rad2deg(in.launch.lat), rad2deg(in.launch.lon), 0];

%%
% % --- UTC 시간축 생성 ---
% t_datetime = in.T0 + seconds(t);
% 
% % --- 태양조명조건 계산 ---
% SAA_th = 10;
% results = STEP3_Solar_visibility(in.target, in.Tmission, SAA_th);
% above_periods = STEP3_findAboveThresholdPeriods(results);

% % --- 조명여부 판정 ---
% isLit = false(size(t_datetime));
% for k = 1:height(above_periods)
%     inPeriod = (t_datetime >= above_periods.Start_UTC(k)) & ...
%                (t_datetime <= above_periods.End_UTC(k));
%     isLit = isLit | inPeriod;
% end
% % --- 가시성과 조명 동시 만족 ---
% vis_final = vis_mission & isLit;
% vis_idx = vis_final(:);
% 
% % 연속 구간 식별
% diff_idx  = diff([0; vis_idx; 0]);
% start_idx = find(diff_idx == 1);
% end_idx   = find(diff_idx == -1) - 1;
% 
% Nseg = numel(start_idx);
% 
% % 1) 각 컬럼을 타임존 포함으로 미리 생성
% Start_UTC    = NaT(Nseg,1,'TimeZone','UTC');
% End_UTC      = NaT(Nseg,1,'TimeZone','UTC');
% Start_Local  = NaT(Nseg,1,'TimeZone','Asia/Seoul');
% End_Local    = NaT(Nseg,1,'TimeZone','Asia/Seoul');
% Duration_sec = zeros(Nseg,1);
% % 2) 채우기
% for kk = 1:Nseg
%     t1 = t_datetime(start_idx(kk));  % t_datetime 자체가 UTC 타임존이면 OK
%     t2 = t_datetime(end_idx(kk));
% 
%     % 로컬 변환
%     t1_local = t1; t1_local.TimeZone = 'Asia/Seoul';
%     t2_local = t2; t2_local.TimeZone = 'Asia/Seoul';
% 
%     Start_UTC(kk)    = t1;
%     End_UTC(kk)      = t2;
%     Start_Local(kk)  = t1_local;
%     End_Local(kk)    = t2_local;
%     Duration_sec(kk) = seconds(t2 - t1);
% end
% 
% % 3) 테이블로 묶기 (이 시점엔 모든 컬럼의 타임존이 이미 설정됨)
% vis_summary = table(Start_UTC, End_UTC, Start_Local, End_Local, Duration_sec);
% 
% disp('==============================================');
% disp('Visible + Illuminated Periods (Satellite Perspective)');
% disp(vis_summary);


%%
% figure('Name','Ground Track','Color','w');
% % set(gcf,'OuterPosition',[100 1200 800 400]);
% latlim = [-90  90];
% lonlim = [  0 360];
% ax = worldmap(latlim,lonlim);
% [N,R] = egm96geoid;
% h = geoshow(ax,N,R,'DisplayType','surface',...
%             'FaceColor',[215 238 237]/255);
% h.ZData = h.ZData - 30;                          % 해수면 ↘︎ 30 m
% land = readgeotable("landareas.shp");
% geoshow(ax,land,"FaceColor",[191 226 202]/255,...
%               "EdgeColor","none");
% 
% 
% hold on
% title(sprintf('Ground Track  |  a = %.1f km   i = %.2f°', a,rad2deg(inc)));
% xlabel('Longitude (deg)');
% ylabel('Latitude (deg)');
% 
% 
% 
% hTrack = geoshow(GroundHist(:,1),GroundHist(:,2), ...
%                  'DisplayType','Line','Color','k');
% 
% hVis = highlight_visible_on_map(ax, GroundHist(:,1:2), vis_mission);
% % hVis2 = highlight_visible_on_map(ax, GroundHist(:,1:2), vis_mission2,1);
% 
% 
% % hVis_lit = highlight_visible_on_map(ax, GroundHist(:,1:2), vis_final, 1);
% 
% 
% 
% 
% hT   = geoshow(tlat,tlon, ...
%                'DisplayType','Point','Marker','s','MarkerEdgeColor','none', ...
%                'MarkerFaceColor',[0.85 0.33 0.10],'MarkerSize',15);   % ■ 빨강
% hTCircle_C1 = geoshow(ax, tlatcC1, tloncC1, 'DisplayType','line', 'Color','k', ...
%                 'LineWidth',1.5);                               % ■ 파란색원
% hTCircle_C2 = geoshow(ax, slatcC2, sloncC2, 'DisplayType','line', 'Color','k', ...
%                 'LineWidth',1.5);                               % ■ 파란색원
% 
% % hTCircle_C3 = geoshow(ax, tlatcC3, tloncC3, 'DisplayType','line', 'Color','b', ...
% %                 'LineWidth',1.5);                               % ■ 파란색원
% hLCircle = geoshow(ax, llatc, llonc, 'DisplayType','line', 'Color','b', ...
%                 'LineWidth',1.5);                               % ■ 파란색원
% 
% 
% hL   = geoshow(rad2deg(in.launch.lat),rad2deg(in.launch.lon), ...
%                'DisplayType','Point','Marker','s','MarkerEdgeColor','none', ...
%                'MarkerFaceColor','b','MarkerSize',15);   % ■ 파랑
% hS   = geoshow(X0lla(1),X0lla(2), ...
%                'DisplayType','Point','Marker','s','MarkerEdgeColor','none', ...
%                'MarkerFaceColor','k','MarkerSize',6);               % ■ 블랙
% hE   = geoshow(Xflla(1),Xflla(2), ...
%                'DisplayType','Point','Marker','s','MarkerEdgeColor','none', ...
%                'MarkerFaceColor','g','MarkerSize',6);               % ■ 초록
% 
% % hCircle = geoshow(ax, latc, lonc, 'DisplayType','line', 'Color','k', ...
% %                 'LineWidth',1.5);                               % ■ 파란색원
% 
% legend([hTrack,  hT,hL, hS,hE], ...
%        {'Ground track','Target','Launch Site','Orbit Injection','End'}, ...
%        'Location','southoutside','Orientation','horizontal');
% 
% % 폰트 일괄 설정
% set(findall(gcf,'-property','FontSize'),...
%     'FontName','Times New Roman','FontSize',12);

%%
% zoomMargin = 3; % deg
% figure('Name','Ground Track (Zoomed near Target)','Color','w');
% 
% % 목표 중심으로 worldmap 생성
% latlim = [tlat-zoomMargin, tlat+zoomMargin];
% lonlim = [tlon-zoomMargin, tlon+zoomMargin];
% ax = worldmap(latlim, lonlim);
% 
% [N,R] = egm96geoid;
% h = geoshow(ax, N, R, 'DisplayType','surface', ...
%             'FaceColor',[215 238 237]/255);
% h.ZData = h.ZData - 30;
% 
% land = readgeotable("landareas.shp");
% geoshow(ax, land, "FaceColor",[191 226 202]/255, "EdgeColor","none");
% 
% hold on
% % 이후 GroundHist, Target, Launch 등 기존 geoshow 코드 그대로
% title(sprintf('Ground Track  |  a = %.1f km   i = %.2f°', a,rad2deg(inc)));
% xlabel('Longitude (deg)');
% ylabel('Latitude (deg)');
% 
% hTrack = geoshow(GroundHist(:,1),GroundHist(:,2), ...
%                  'DisplayType','Line','Color','k');
% 
% hVis = highlight_visible_on_map(ax, GroundHist(:,1:2), vis_mission);
% % hVis2 = highlight_visible_on_map(ax, GroundHist(:,1:2), vis_mission2,1);
% 
% 
% hT   = geoshow(tlat,tlon, ...
%                'DisplayType','Point','Marker','s','MarkerEdgeColor','none', ...
%                'MarkerFaceColor',[0.85 0.33 0.10],'MarkerSize',15);   % ■ 빨강
% hTCircle_C1 = geoshow(ax, tlatcC1, tloncC1, 'DisplayType','line', 'Color','k', ...
%                 'LineWidth',1.5);                               % ■ 파란색원
% % hTCircle_C2 = geoshow(ax, slatcC2, sloncC2, 'DisplayType','line', 'Color','k', ...
% %                 'LineWidth',1.5);                               % ■ 파란색원
% 
% % hTCircle_C3 = geoshow(ax, tlatcC3, tloncC3, 'DisplayType','line', 'Color','b', ...
% %                 'LineWidth',1.5);                               % ■ 파란색원
% % hLCircle = geoshow(ax, llatc, llonc, 'DisplayType','line', 'Color','b', ...
% %                 'LineWidth',1.5);                               % ■ 파란색원
% 
% 
% hL   = geoshow(rad2deg(in.launch.lat),rad2deg(in.launch.lon), ...
%                'DisplayType','Point','Marker','s','MarkerEdgeColor','none', ...
%                'MarkerFaceColor','b','MarkerSize',15);   % ■ 파랑
% hS   = geoshow(X0lla(1),X0lla(2), ...
%                'DisplayType','Point','Marker','s','MarkerEdgeColor','none', ...
%                'MarkerFaceColor','k','MarkerSize',6);               % ■ 블랙
% hE   = geoshow(Xflla(1),Xflla(2), ...
%                'DisplayType','Point','Marker','s','MarkerEdgeColor','none', ...
%                'MarkerFaceColor','g','MarkerSize',6);               % ■ 초록
% 
% % hCircle = geoshow(ax, latc, lonc, 'DisplayType','line', 'Color','k', ...
% %                 'LineWidth',1.5);                               % ■ 파란색원
% 
% 
% % 폰트 일괄 설정
% set(findall(gcf,'-property','FontSize'),...
%     'FontName','Times New Roman','FontSize',12);

%%
zoomMargin = 3; % deg
figure('Name','Ground Track (Zoomed near Target)','Color','w');

% 목표 중심으로 worldmap 생성
latlim = [tlat-zoomMargin, tlat+zoomMargin];
lonlim = [tlon-zoomMargin, tlon+zoomMargin];
ax = worldmap(latlim, lonlim);

[N,R] = egm96geoid;
h = geoshow(ax, N, R, 'DisplayType','surface', ...
            'FaceColor',[215 238 237]/255);
h.ZData = h.ZData - 30;

land = readgeotable("landareas.shp");
geoshow(ax, land, "FaceColor",[191 226 202]/255, "EdgeColor","none");

hold on
% 이후 GroundHist, Target, Launch 등 기존 geoshow 코드 그대로
title(sprintf('Ground Track  |  a = %.1f km   i = %.2f°', a,rad2deg(inc)));
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

hTrack = geoshow(GroundHist(:,1),GroundHist(:,2), ...
                 'DisplayType','Line','Color','k');

% hVis = highlight_visible_on_map(ax, GroundHist(:,1:2), vis_mission);
% hVis2 = highlight_visible_on_map(ax, GroundHist(:,1:2), vis_mission2,1);
% hVis_lit = highlight_visible_on_map(ax, GroundHist(:,1:2), vis_final, 1);


hT   = geoshow(tlat,tlon, ...
               'DisplayType','Point','Marker','s','MarkerEdgeColor','none', ...
               'MarkerFaceColor',[0.85 0.33 0.10],'MarkerSize',15);   % ■ 빨강
hTCircle_C1 = geoshow(ax, tlatcC1, tloncC1, 'DisplayType','line', 'Color','k', ...
                'LineWidth',1.5);                               % ■ 파란색원
% hTCircle_C2 = geoshow(ax, slatcC2, sloncC2, 'DisplayType','line', 'Color','k', ...
%                 'LineWidth',1.5);                               % ■ 파란색원

% hTCircle_C3 = geoshow(ax, tlatcC3, tloncC3, 'DisplayType','line', 'Color','b', ...
%                 'LineWidth',1.5);                               % ■ 파란색원
% hLCircle = geoshow(ax, llatc, llonc, 'DisplayType','line', 'Color','b', ...
%                 'LineWidth',1.5);                               % ■ 파란색원


hL   = geoshow(rad2deg(in.launch.lat),rad2deg(in.launch.lon), ...
               'DisplayType','Point','Marker','s','MarkerEdgeColor','none', ...
               'MarkerFaceColor','b','MarkerSize',15);   % ■ 파랑
hS   = geoshow(X0lla(1),X0lla(2), ...
               'DisplayType','Point','Marker','s','MarkerEdgeColor','none', ...
               'MarkerFaceColor','k','MarkerSize',6);               % ■ 블랙
hE   = geoshow(Xflla(1),Xflla(2), ...
               'DisplayType','Point','Marker','s','MarkerEdgeColor','none', ...
               'MarkerFaceColor','g','MarkerSize',6);               % ■ 초록

% hCircle = geoshow(ax, latc, lonc, 'DisplayType','line', 'Color','k', ...
%                 'LineWidth',1.5);                               % ■ 파란색원


% 폰트 일괄 설정
set(findall(gcf,'-property','FontSize'),...
    'FontName','Times New Roman','FontSize',12);

%%
% % Plot orbit in ECEF reference
% figure('Name','Orbit ECEF','Color','w');
% 
% % --- 3D 지구 그리기 (earth_sphere 함수 호출) ---
% earth_sphere('m');
% hold on;
% 
% % 궤도 플롯
% plot3(Eph_ecef(:,2), Eph_ecef(:,3), Eph_ecef(:,4), 'k', 'LineWidth', 1);
% 
% % --- 시작/종료/목표/발사 지점 좌표 계산 (ECEF, m) ---
% % (이 부분은 이전 코드와 동일합니다)
% start_ecef = Eph_ecef(1, 2:4);
% end_ecef   = Eph_ecef(end, 2:4);
% target_ecef_m = lla2ecef([tlat, tlon, target_lladeg(3)*100]);
% launch_ecef_m = lla2ecef([rad2deg(in.launch.lat), rad2deg(in.launch.lon), in.launch.alt]);
% 
% 
% 
% 
% % --- 지점 표시 ---
% hT_3d = plot3(target_ecef_m(1), target_ecef_m(2), target_ecef_m(3), ...
%               's','MarkerFaceColor',[0.85 0.33 0.10],'MarkerSize',15,'MarkerEdgeColor','none'); % 목표(T)
% hL_3d = plot3(launch_ecef_m(1), launch_ecef_m(2), launch_ecef_m(3), ...
%               's','MarkerFaceColor','y','MarkerSize',15,'MarkerEdgeColor','none'); % 발사(L)
% hS_3d = plot3(start_ecef(1), start_ecef(2), start_ecef(3), ...
%               's','MarkerFaceColor','k','MarkerSize',8,'MarkerEdgeColor','none'); % 전파 시작(S)
% hE_3d = plot3(end_ecef(1), end_ecef(2), end_ecef(3), ...
%               's','MarkerFaceColor','g','MarkerSize',8,'MarkerEdgeColor','none'); % 전파 종료(E)
% scale = 1.2;  % 화살표 길이 배율 (지구 반지름보다 조금 크게)
% quiver3(0,0,0, ...
%         target_ecef_m(1)*scale, ...
%         target_ecef_m(2)*scale, ...
%         target_ecef_m(3)*scale, ...
%         0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
% 
% 
% 
% % 그래프 설정
% grid on;
% axis equal;
% title('Orbit in ECEF (Earth-Fixed) Frame (m)');
% xlabel('X (m)');
% ylabel('Y (m)');
% zlabel('Z (m)');
% legend([hT_3d, hL_3d, hS_3d, hE_3d], 'Target', 'Launch', 'Start', 'End', 'Location', 'best');
% % campos(target_ecef_m*100);      % 카메라 위치를 타겟의 ECEF 좌표로 이동
% % camtarget([0 0 0]);         % 카메라가 바라볼 지점을 지구 중심(원점)으로 설정
% % camup([0 0 1]);             % 카메라의 위쪽 방향을 Z축으로 설정 (지구 북극 방향)
% % 
% %%
% % Plot orbit in ECI reference
% figure('Name','Orbit ECI','Color','w');
% 
% % --- 3D 지구 그리기 (earth_sphere 함수 호출) ---
% % 단위 'm'를 지정하여 미터 단위의 지구를 그립니다.
% earth_sphere('m'); 
% hold on;
% 
% % 궤도 플롯
% plot3(Eph_eci(:,2), Eph_eci(:,3), Eph_eci(:,4), 'r', 'LineWidth', 2);
% 
% % 그래프 설정
% grid on;
% axis equal;
% title('Orbit in ECI (Inertial) Frame (m)');
% xlabel('X (m)');
% ylabel('Y (m)');
% zlabel('Z (m)');
% view(3);


%%
%  a, ecc, incl, omega, argp, nu, m
% figure;
% subplot(4,2,1)
% plot(Eph_osc(:,3))
% subplot(4,2,2)
% plot(Eph_osc(:,4))
% subplot(4,2,3)
% plot(Eph_osc(:,5))
% subplot(4,2,4)
% plot(Eph_osc(:,6))
% subplot(4,2,5)
% plot(Eph_osc(:,7))
% subplot(4,2,6)
% plot(Eph_osc(:,8))
% subplot(4,2,7)
% plot(Eph_osc(:,9))

% %%
% time = Eph_eci(:,1);
% r_mag = sqrt(Eph_eci(:,2).^2 + Eph_eci(:,3).^2 + Eph_eci(:,4).^2);
% p = polyfit(time, r_mag, 2); % 2차 추세선
% r_trend = polyval(p, time);
% 
% figure;
% plot(time/3600, r_mag/1e3, 'Color',[0.7 0.7 0.7]); % 원 데이터
% hold on;
% plot(time/3600, r_trend/1e3, 'r-', 'LineWidth',1.5); % 추세선
% grid on;
% xlabel('Time [hours]');
% ylabel('r magnitude [km]');
% title('Geocentric distance with polynomial trend');
% legend('Raw r(t)','Polynomial trend');
% 











function [eta_deg, vis_horizon, vis_mission] = fov_nadir_conical( ...
    r_ecef_km, target_lladeg, Re_km, eta_thresh_deg, useAccurateHorizon)

% FOV (look angle) for a nadir-pointing, axisymmetric conical sensor.
% Inputs
%   r_ecef_km         : [N x 3] or [3 x N] satellite ECEF positions [km]
%   target_lladeg     : [lat_deg, lon_deg, alt_km] (alt_km=0 for ground)
%   Re_km             : Earth radius [km] (e.g., 6378.137)
%   eta_thresh_deg    : mission threshold (e.g., 10)
%   useAccurateHorizon: true = use eta_h exact; false = coarse check
% Outputs
%   eta_deg     : [N x 1] look angle history [deg]
%   vis_horizon : [N x 1] horizon visibility (true if not occulted)
%   vis_mission : [N x 1] (eta < eta_thresh) & vis_horizon
%
% -----------이론
%
% 지구중심 - 위성 - 목표지점 삼각형을 생각 (O-S-T)
% O의 각: 지상 거리(Λ)
% S의 각: 조준각(Boresight Angle, η)
% T의 각: 고도각(elevation Angle, el) +90deg 
% 경사 거리(Slant Range, ρ):  위성에서 목표 지점까지의 실제 직선거리
% FOV : 센서가 한 번에 관측할 수 있는 총 영역
    arguments
        r_ecef_km double
        target_lladeg (1,3) double
        Re_km (1,1) double {mustBePositive}
        eta_thresh_deg (1,1) double {mustBeNonnegative} = 10
        useAccurateHorizon (1,1) logical = true
    end

    % Ensure shape = [N x 3]
    if size(r_ecef_km,2) ~= 3
        r_ecef_km = r_ecef_km.'; % assume [3 x N]
    end
    N = size(r_ecef_km,1);

    % Target ECEF (sphere) in km
    lat = target_lladeg(1);
    lon = target_lladeg(2);
    alt = target_lladeg(3);
    rTP_ecef_km = lla2ecef([lat,lon,alt]);           % [1 x 3] [m]
    rTP_ecef_km = rTP_ecef_km./1e3;

    % Unit vectors
    rs_km  = sqrt(sum(r_ecef_km.^2,2));            % ||r_sat||
    rhat   = r_ecef_km ./ rs_km;                    % r̂_sat
    rTPhat = rTP_ecef_km ./ norm(rTP_ecef_km);      % r̂_TP

    % Ground-range central angle Λ
    % 위성 나디르방향 단위벡터와 타겟방향 단위벡터 내적 
    cosL = -rhat * rTPhat.';                        % [N x 1]
    cosL = max(-1,min(1,cosL));                     % clamp
    sinL = sqrt(max(0, 1 - cosL.^2));

    % Slant range ρ (target at surface + alt)
    rho_km = sqrt(sum( (r_ecef_km - rTPhat.*(Re_km+(alt/1e3))).^2, 2 ));

    % Look angle η (atan2 form, numerically robust)
    
    num = Re_km .* sinL;                            % R*sinΛ
    den = rs_km - Re_km .* cosL;                    % r_s - R*cosΛ
    eta_rad = atan2(num, den);
    eta_deg = rad2deg(eta_rad);

    % Horizon check
    if useAccurateHorizon
        % eta_h = asin(R/rs)
        ratio = max(0, min(1, Re_km ./ rs_km));
        eta_h_rad = asin(ratio);
        vis_horizon = eta_rad <= eta_h_rad + 1e-12;
    else
        % coarse: rho < R (loose and not exact)
        vis_horizon = rho_km < Re_km;
    end

    % Mission threshold
    vis_mission = (eta_deg < eta_thresh_deg) & vis_horizon;
end


function hVis = highlight_visible_on_map(ax, GroundHist_latlon, vis_mission, varargin)
% GroundHist_latlon: [N x 2] [lat_deg, lon_deg]
% vis_mission      : [N x 1] logical

    if nargin < 3 || isempty(vis_mission)
        vis_mission = true(size(GroundHist_latlon,1),1);
    end

    vis_mission = logical(vis_mission(:));

    lat = GroundHist_latlon(:,1);
    lon = GroundHist_latlon(:,2);

    % 1) 비가시 구간을 NaN으로 만들어 선을 끊는다
    lat(~vis_mission) = NaN;
    lon(~vis_mission) = NaN;

    % 2) 경도 점프(지도 래핑)에서도 끊기 (옵션)
    %    wrapTo180 사용 시 변화가 큰 지점에서 끊어주면 긴 대각선 방지
    lon180 = wrapTo180(lon);
    jump = [false; abs(diff(lon180)) > 180];   % 혹은  > 120~150 정도로 느슨하게
    lat(jump) = NaN;
    lon(jump) = NaN;

    hold(ax,'on');
    if ~isempty(varargin) && varargin{1} == 1
        hVis = geoshow(ax, lat, lon, ...
        'DisplayType','line', 'Color','b', 'LineWidth',3);

    else
        hVis = geoshow(ax, lat, lon, ...
        'DisplayType','line', 'Color','r', 'LineWidth',3);
    end
end
