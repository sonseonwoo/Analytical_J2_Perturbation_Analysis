%% ─── Mean (sol) 궤도 적분 및 지상궤적 생성 ─────────────────────────
a          = best.sol.a*1e-3;
inc        = best.sol.i;
RAAN       = best.sol.RAAN0;

if isfield(best.sol, 'u0')
    arglat     = best.sol.u0;
    [r0_ijk,v0_ijk]  = coe2rv(a,0,inc,RAAN,0,0,arglat,0,0); 
else
    e   = best.sol.e;
    AoP = best.sol.AoP0;
    nu  = best.sol.nu0;
    [r0_ijk,v0_ijk]  = coe2rv(a,e,inc,RAAN,AoP,nu,0,0,0); 
end

% Propagation setting
in.n       = 70;
in.m       = 70;
in.sun     = 1; in.moon = 1; in.planets = 1;
in.sRad    = 1; in.drag  = 1;
in.SolidEarthTides = 1; in.OceanTides = 1; in.Relativity = 1;

% Initial state [m]
Y0 = [r0_ijk; v0_ijk] * 1e3;
Step   = 1;   
N_Step = (best.tp);
[t, yout] = ode113(@Accel, (0:Step:N_Step*Step), Y0, in.options);

% Process orbit data
[Eph_eci, Eph_ecef, Eph_osc_dummy] = processOrbitData(yout, t);

% % Ground track (mean)
% earthSphere = referenceSphere('Earth','meter');
% x = Eph_ecef(:,2); y = Eph_ecef(:,3); z = Eph_ecef(:,4);
% [lat, lon, alt] = ecef2geodetic(x, y, z, earthSphere);
% GroundHist = [rad2deg(lat), mod(rad2deg(lon),360), alt];
% 
GroundHist = ecef2lla(Eph_ecef(:,2:4));


% ─── Osc 궤도 적분 및 지상궤적 생성 ─────────────────────────
a_osc    = best.osc.a*1e-3;
inc_osc  = best.osc.i;
RAAN_osc = best.osc.RAAN;

if isfield(best.osc, 'u')
    arglat = best.osc.u;
    [r0_ijk,v0_ijk]  = coe2rv(a_osc,0,inc_osc,RAAN_osc,0,0,arglat,0,0); 
else
    e   = best.osc.e;
    AoP = best.osc.AoP0;
    nu  = best.osc.nu0;
    [r0_ijk,v0_ijk]  = coe2rv(a_osc,e,inc_osc,RAAN_osc,AoP,nu,0,0,0); 
end

Y0_osc = [r0_ijk; v0_ijk] * 1e3;
[t_osc, yout_osc] = ode113(@Accel, (0:Step:N_Step*Step), Y0_osc, in.options);

% Process orbit data for osc
[Eph_eci_osc, Eph_ecef_osc, Eph_osc] = processOrbitData(yout_osc, t_osc);

% % Ground track (osc)
% x = Eph_ecef_osc(:,2); y = Eph_ecef_osc(:,3); z = Eph_ecef_osc(:,4);
% [lat, lon, alt] = ecef2geodetic(x, y, z, earthSphere);
% GroundHist_osc = [rad2deg(lat), mod(rad2deg(lon),360), alt];

GroundHist_osc = ecef2lla(Eph_ecef_osc(:,2:4));


%%
error_oscVSmean(Eph_ecef, Eph_ecef_osc, GroundHist, GroundHist_osc, in)
%%

