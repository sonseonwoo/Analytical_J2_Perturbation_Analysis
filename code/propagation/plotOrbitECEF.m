function plotOrbitECEF(Eph_ecef, target_lla, launch_lla)
% ECEF(지구고정) 좌표계 기준 3D 위성 궤도와 애니메이션을 그립니다.
% Inputs:
%   Eph_ecef   : [Nx8] ECEF 데이터 [t, r, v, alt] (m)
%   C          : 상수 구조체 (RE 포함)
%   target_lla : [1x3] 타겟 위도, 경도(deg), 고도(m)
%   launch_lla : [1x3] 발사지 위도, 경도(deg), 고도(m)

figure('Name','Orbit ECEF','Color','w');
earth_sphere('m');
hold on;

% 전체 궤적 미리보기
plot3(Eph_ecef(:,2), Eph_ecef(:,3), Eph_ecef(:,4), 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
% plot3(Eph_ecef(:,2), Eph_ecef(:,3), Eph_ecef(:,4), 'Color', 'r', 'LineWidth', 10);


% 주요 지점 계산 및 표시
start_ecef = Eph_ecef(1, 2:4);
end_ecef   = Eph_ecef(end, 2:4);
target_ecef_m = lla2ecef(target_lla);
launch_ecef_m = lla2ecef(launch_lla);

hT = plot3(target_ecef_m(1), target_ecef_m(2), target_ecef_m(3), 's','MarkerFaceColor',[0.85 0.33 0.10],'MarkerSize',15,'MarkerEdgeColor','none');
hL = plot3(launch_ecef_m(1), launch_ecef_m(2), launch_ecef_m(3), 's','MarkerFaceColor','b','MarkerSize',15,'MarkerEdgeColor','none');
hS = plot3(start_ecef(1), start_ecef(2), start_ecef(3), 's','MarkerFaceColor','k','MarkerSize',8,'MarkerEdgeColor','none');
hE = plot3(end_ecef(1), end_ecef(2), end_ecef(3), 's','MarkerFaceColor','g','MarkerSize',8,'MarkerEdgeColor','none');

% 그래프 설정
grid on; axis equal;
title('Orbit in ECEF (Earth-Fixed) Frame (m)');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
legend([hT, hL, hS, hE], 'Target', 'Launch', 'Start', 'End', 'Location', 'best');
% view(3);
campos(target_ecef_m);      % 카메라 위치를 타겟의 ECEF 좌표로 이동
camtarget([0 0 0]);         % 카메라가 바라볼 지점을 지구 중심(원점)으로 설정
camup([0 0 1]);             % 카메라의 위쪽 방향을 Z축으로 설정 (지구 북극 방향)


% 애니메이션
an_3d_track = animatedline(gca, 'LineStyle','-', 'Color','k', 'LineWidth',2);
an_3d_sat = animatedline(gca, 'Marker','o', 'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'MarkerSize',10);

for k = 1:size(Eph_ecef, 1)
    addpoints(an_3d_track, Eph_ecef(k,2), Eph_ecef(k,3), Eph_ecef(k,4));
    clearpoints(an_3d_sat);
    addpoints(an_3d_sat, Eph_ecef(k,2), Eph_ecef(k,3), Eph_ecef(k,4));
    % drawnow limitrate;
    pause(0.01)
end

end