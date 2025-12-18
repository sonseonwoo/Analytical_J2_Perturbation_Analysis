function err = error_oscVSmean(Eph_ecef, Eph_ecef_osc, GroundHist, GroundHist_osc, in)
% ===============================================================
% error_oscVSmean
%   Mean(=sol) / Osc 궤도의 ECEF 및 GroundHist 기반 오차 계산
%
% INPUTS
%   Eph_ecef        : mean 궤도의 ECEF 위치 벡터 [N×4] (t, x, y, z)
%   Eph_ecef_osc    : osc 궤도의 ECEF 위치 벡터 [N×4]
%   GroundHist      : mean 궤도의 [lat_deg, lon_deg, alt_m]
%   GroundHist_osc  : osc 궤도의 [lat_deg, lon_deg, alt_m]
%   in              : 구조체 (in.target.lat/lon 포함)
%
% OUTPUTS
%   err 구조체 (모든 값은 [m])
% ===============================================================

    % === 1. 타겟 좌표 (ECEF, [m]) ===
    tgt_lat  = rad2deg(in.target.lat);
    tgt_lon  = rad2deg(in.target.lon);
    tgt_alt  = 0;
    r_target = lla2ecef([tgt_lat, tgt_lon, tgt_alt]);  % [m]

    % === 2. Mean / Osc 위치 벡터 ===
    r_start_mean = Eph_ecef(1,  2:4);
    r_final_mean = Eph_ecef(end,2:4);
    r_start_osc  = Eph_ecef_osc(1,  2:4);
    r_final_osc  = Eph_ecef_osc(end,2:4);

    % === 3. ECEF 거리 기반 norm 비교 ===
    err.ecef.mean_vs_osc_start = norm(r_start_mean - r_start_osc);
    err.ecef.mean_vs_osc_end   = norm(r_final_mean - r_final_osc);
    err.ecef.mean_vs_target_end = norm(r_final_mean - r_target);
    err.ecef.osc_vs_target_end  = norm(r_final_osc - r_target);

    % === 4. GroundHist 기반 거리 비교 ===
    % (1) Mean–Osc start / end
    arc_deg_start = distance(GroundHist(1,1), GroundHist(1,2), ...
                             GroundHist_osc(1,1), GroundHist_osc(1,2));
    arc_deg_end   = distance(GroundHist(end,1), GroundHist(end,2), ...
                             GroundHist_osc(end,1), GroundHist_osc(end,2));
    err.ground.mean_vs_osc_start = deg2km(arc_deg_start) * 1e3;  % [m]
    err.ground.mean_vs_osc_end   = deg2km(arc_deg_end)   * 1e3;  % [m]

    % (2) Target–Mean / Target–Osc
    tgt_lat = rad2deg(in.target.lat);
    tgt_lon = rad2deg(in.target.lon);
    arc_deg_target_mean = distance(tgt_lat, tgt_lon, ...
                                   GroundHist(end,1), GroundHist(end,2));
    arc_deg_target_osc = distance(tgt_lat, tgt_lon, ...
                                  GroundHist_osc(end,1), GroundHist_osc(end,2));
    err.ground.mean_vs_target_end = deg2km(arc_deg_target_mean) * 1e3;  % [m]
    err.ground.osc_vs_target_end  = deg2km(arc_deg_target_osc)  * 1e3;  % [m]

    % === 5. 출력 ===
    fprintf('\n[ ECEF Position Errors ]\n');
    fprintf('  Mean–Osc  start : %.6f m\n', err.ecef.mean_vs_osc_start);
    fprintf('  Mean–Osc  end   : %.6f m\n', err.ecef.mean_vs_osc_end);
    fprintf('  Mean–Target end : %.6f m\n', err.ecef.mean_vs_target_end);
    fprintf('  Osc–Target  end : %.6f m\n', err.ecef.osc_vs_target_end);

    fprintf('\n[ GroundTrack Distance Errors ]\n');
    fprintf('  Mean–Osc  start : %.6f m\n', err.ground.mean_vs_osc_start);
    fprintf('  Mean–Osc  end   : %.6f m\n', err.ground.mean_vs_osc_end);
    fprintf('  Mean–Target end : %.6f m\n', err.ground.mean_vs_target_end);
    fprintf('  Osc–Target  end : %.6f m\n', err.ground.osc_vs_target_end);

end
