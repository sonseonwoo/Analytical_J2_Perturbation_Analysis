function saveEphemerisFiles(Eph_eci, Eph_ecef, Eph_osc)
    % processOrbitData 함수로 계산된 결과물들을 3개의 텍스트 파일로 저장합니다.
    global in
    Mjd_UTC = in.Mjd_UTC;
    n_steps = size(Eph_eci, 1);

    % --- 1. ECI 데이터 저장 ---
    fid = fopen('SatelliteStates_eci.txt','w');
    fprintf(fid, '# ECI State Vectors (Units: m, m/s)\n');
    fprintf(fid, '%-25s %14s %14s %18s %14s %14s %14s %14s\n', ...
        '# Epoch (UTCG)', 'x(m)', 'y(m)', 'z(m)', 'vx(m/s)', 'vy(m/s)', 'vz(m/s)', 'a(m)');
    for i = 1:n_steps
        [year,month,day,hour,minute,sec] = invjday(Mjd_UTC + Eph_eci(i,1)/86400);
        fprintf(fid,'  %4d/%02d/%02d  %02d:%02d:%06.3f',year,month,day,hour,minute,sec);
        fprintf(fid,'  %16.3f%16.3f%16.3f%14.3f%14.3f%14.3f%16.3f\n', Eph_eci(i,2:end));
    end
    fclose(fid);

    % --- 2. ECEF 데이터 저장 ---
    fid = fopen('SatelliteStates_ecef.txt','w');
    fprintf(fid, '# ECEF State Vectors (Units: m, m/s)\n');
    fprintf(fid, '%-25s %14s %14s %18s %14s %14s %14s %14s\n', ...
        '# Epoch (UTCG)', 'x(m)', 'y(m)', 'z(m)', 'vx(m/s)', 'vy(m/s)', 'vz(m/s)', 'alt(m)');
    for i = 1:n_steps
        [year,month,day,hour,minute,sec] = invjday(Mjd_UTC + Eph_ecef(i,1)/86400);
        fprintf(fid,'  %4d/%02d/%02d  %02d:%02d:%06.3f',year,month,day,hour,minute,sec);
        fprintf(fid,'  %16.3f%16.3f%16.3f%14.3f%14.3f%14.3f%16.3f\n', Eph_ecef(i,2:end));
    end
    fclose(fid);

    % --- 3. 궤도 요소 데이터 저장 ---
    fid = fopen('SatelliteState_coe.txt', 'w');
    fprintf(fid, '# Satellite Osculating Orbital Elements (Units: Angles in Degrees, Distances in km)\n');
    fprintf(fid, '%-25s %14s %14s %14s %12s %12s %12s %12s %12s %12s %12s %12s %14s %14s\n', ...
    '# Epoch (UTCG)', 'p(km)', 'a(km)', 'ecc', 'incl(deg)', 'RAAN(deg)', 'argp(deg)', 'nu(deg)', 'M(deg)', 'arglat(deg)', 'truelon(deg)', 'lonper(deg)', 'xi', 'xi_J2');
    for i = 1:n_steps
        [year,month,day,hour,minute,sec] = invjday(Mjd_UTC + Eph_osc(i,1)/86400);
        time_str = sprintf('%4d/%02d/%02d %02d:%02d:%06.3f', year,month,day,hour,minute,sec);
        data_row = Eph_osc(i, 2:end);
        angle_indices = 4:11;
        for k = angle_indices
            if ~isnan(data_row(k))
                data_row(k) = rad2deg(data_row(k));
            end
        end
        fprintf(fid, '%-25s %14.4f %14.4f %14.9f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %14.6f %14.6f\n', time_str, data_row);
    end
    fclose(fid);
end