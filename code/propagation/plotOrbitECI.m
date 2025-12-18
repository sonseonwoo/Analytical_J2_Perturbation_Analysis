function plotOrbitECI(Eph_eci)
% ECI(관성) 좌표계 기준 3D 위성 궤도를 그립니다.
% Inputs:
%   Eph_eci : [Nx8] ECI 데이터 [t, r, v, a] (m)
%   C       : 상수 구조체 (RE 포함)
figure('Name','Orbit ECI','Color','w');

earth_sphere('m'); % 단위 'm'
hold on;

plot3(Eph_eci(:,2), Eph_eci(:,3), Eph_eci(:,4), 'r', 'LineWidth', 2);

grid on; axis equal;
title('Orbit in ECI (Inertial) Frame (m)');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
view(3);

end