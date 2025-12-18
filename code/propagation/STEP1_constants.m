% ================================
% constants.m
% ================================
function C = STEP1_constants()
% Physical & conversion constants for FAO–Revisit mixed design
%
% OUTPUT
%   C : struct with fields below (all SI base units)


% ───  E arth & gravity  ───────────────────────────────────────────
C.muE   = 3.986004418e14;            % [m^3/s^2]  gravitational parameter μ_e
C.RE    = 6378136.3;                 % [m]        mean equatorial radius R_e
C.J2    = 1.08262668e-3;             % [‒]        2nd zonal harmonic coefficient
C.f     = 1/298.257223563;           % [-] 편평도 


% ───  Earth rotation  ─────────────────────────────────────────────
C.wE    = 7.2921150e-5;              % [rad/s]    sidereal rotation rate (≈15.041°/h)

% From Mahooti

% Mathematical constants 
C.pi2       = 2*pi;                % 2pi
C.Rad       = pi/180;              % Radians per degree
C.Deg       = 180/pi;              % Degrees per radian 
C.Arcs      = 3600*180/pi;         % Arcseconds per radian

% General
C.MJD_J2000 = 51544.5;             % Modified Julian Date of J2000
C.T_B1950   = -0.500002108;        % Epoch B1950
C.c_light   = 299792457.999999984; % Speed of light  [m/s]; DE440
C.AU        = 149597870699.999988; % Astronomical unit [m]; DE440

% Physical parameters of the Earth, Sun and Moon

% Equatorial radius and flattening
C.R_Earth   = 6378.137e3;                          % Earth's radius [m]; WGS-84
C.f_Earth   = 1/298.257223563;                     % Flattening; WGS-84
C.R_Sun     = 696000.0e3;                          % Sun's radius [m]; DE440
C.R_Moon    = 1738.0e3;                            % Moon's radius [m]; DE440

% Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
C.omega_Earth = 15.04106717866910/3600*C.Rad;      % [rad/s]; WGS-84

% Gravitational coefficients
C.GM_Earth   = 398600.4415e9;         			   % [m^3/s^2]; GGM03C & GGM03S
C.GM_Sun     = 132712440041.279419e9; 			   % [m^3/s^2]; DE440
C.GM_Moon    = C.GM_Earth/81.3005682214972154;     % [m^3/s^2]; DE440
C.GM_Mercury = 22031.868551e9; 		  			   % [m^3/s^2]; DE440
C.GM_Venus   = 324858.592000e9;       			   % [m^3/s^2]; DE440
C.GM_Mars    = 42828.375816e9;	      			   % [m^3/s^2]; DE440
C.GM_Jupiter = 126712764.100000e9;    			   % [m^3/s^2]; DE440
C.GM_Saturn  = 37940584.841800e9;     			   % [m^3/s^2]; DE440
C.GM_Uranus  = 5794556.400000e9;      			   % [m^3/s^2]; DE440
C.GM_Neptune = 6836527.100580e9;      			   % [m^3/s^2]; DE440
C.GM_Pluto   = 975.500000e9;	      			   % [m^3/s^2]; DE440



% Solar radiation pressure at 1 AU
C.P_Sol = 1367/C.c_light;                          % [N/m^2] (1367 W/m^2); IERS

% Pi
C.DPI = 3.141592653589793238462643;
% 2Pi
C.D2PI = 6.283185307179586476925287;
% Radians to degrees
C.DR2D = 57.29577951308232087679815;
% Degrees to radians
C.DD2R = 1.745329251994329576923691e-2;
% Radians to arcseconds
C.DR2AS = 206264.8062470963551564734;
% Arcseconds to radians
C.DAS2R = 4.848136811095359935899141e-6;
% Seconds of time to radians
C.DS2R =7.272205216643039903848712e-5;
% Arcseconds in a full circle
C.TURNAS = 1296000.0;
% Milliarcseconds to radians
C.DMAS2R = C.DAS2R / 1e3;
% Length of tropical year B1900 (days)
C.DTY = 365.242198781;
% Seconds per day.
C.DAYSEC = 86400.0;
% Days per Julian year
C.DJY = 365.25;
% Days per Julian century
C.DJC = 36525.0;
% Days per Julian millennium
C.DJM = 365250.0;
% Reference epoch (J2000.0), Julian Date
C.DJ00 = 2451545.0;
% Julian Date of Modified Julian Date zero
C.DJM0 = 2400000.5;
% Reference epoch (J2000.0), Modified Julian Date
C.DJM00 = 51544.5;
% 1977 Jan 1.0 as MJD
C.DJM77 = 43144.0;
% TT minus TAI (s)
C.TTMTAI = 32.184;
% Astronomical unit (m, IAU 2012)
C.DAU = 149597870.7e3;
% Speed of light (m/s)
C.CMPS = 299792458.0;
% Light time for 1 au (s)
C.AULT = C.DAU/C.CMPS;
% Speed of light (au per day)
C.DC = C.DAYSEC/C.AULT;
% L_G = 1 - d(TT)/d(TCG)
C.ELG = 6.969290134e-10;
% L_B = 1 - d(TDB)/d(TCB), and TDB (s) at TAI 1977/1/1.0
C.ELB = 1.550519768e-8;
C.TDB0 = -6.55e-5;
% Schwarzschild radius of the Sun (au)
% = 2 * 1.32712440041e20 / (2.99792458e8)^2 / 1.49597870700e11
C.SRS = 1.97412574336e-8;



end
