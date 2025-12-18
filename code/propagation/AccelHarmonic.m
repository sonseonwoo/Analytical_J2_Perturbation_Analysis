% --------------------------------------------------------------------------
% 
% AccelHarmonic: Computes the acceleration due to the harmonic gravity
%                field of the central body
% 
% Inputs:
%   r           Satellite position vector in the inertial system
%   E           Transformation matrix to body-fixed system
%   Cnm, Snm    Spherical harmonics coefficients (normalized)
%   n_max       Maximum degree
%   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
% 
% Output:
%   a           Acceleration (a=d^2r/dt^2)
% 
% Last modified:   2022/06/16   Meysam Mahooti
% 
% --------------------------------------------------------------------------
% function a = AccelHarmonic(r, E, n_max, m_max)
% 
% global Cnm Snm
% 
% gm    = 398600.4415e9; % [m^3/s^2]; GGM03C & GGM03S
% r_ref = 6378.1363e3;   % Earth's radius [m]; GGM03C & GGM03S
% 
% % Body-fixed position 
% r_bf = E * r;
% 
% % Auxiliary quantities
% d = norm(r_bf);                     % distance
% latgc = asin(r_bf(3)/d);
% lon = atan2(r_bf(2),r_bf(1));
% 
% [pnm, dpnm] = Legendre(n_max,m_max,latgc);
% 
% dUdr = 0;
% dUdlatgc = 0;
% dUdlon = 0;
% q3 = 0; q2 = q3; q1 = q2;
% for n=0:n_max
%     b1 = (-gm/d^2)*(r_ref/d)^n*(n+1);
%     b2 =  (gm/d)*(r_ref/d)^n;
%     b3 =  (gm/d)*(r_ref/d)^n;
%     for m=0:m_max
%         q1 = q1 + pnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
%         q2 = q2 + dpnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
%         q3 = q3 + m*pnm(n+1,m+1)*(Snm(n+1,m+1)*cos(m*lon)-Cnm(n+1,m+1)*sin(m*lon));
%     end
%     dUdr     = dUdr     + q1*b1;
%     dUdlatgc = dUdlatgc + q2*b2;
%     dUdlon   = dUdlon   + q3*b3;
%     q3 = 0; q2 = q3; q1 = q2;
% end
% 
% % Body-fixed acceleration
% r2xy = r_bf(1)^2+r_bf(2)^2;
% 
% ax = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(1)-(1/r2xy*dUdlon)*r_bf(2);
% ay = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1/r2xy*dUdlon)*r_bf(1);
% az =  1/d*dUdr*r_bf(3)+sqrt(r2xy)/d^2*dUdlatgc;
% 
% a_bf = [ax ay az]';
% 
% % Inertial acceleration 
% a = E'*a_bf;

%--------------------------------------------------------------------------
%
% AccelHarmonic: Computes the acceleration due to the harmonic gravity
%                field of the central body. Includes special cases for
%                Two-Body (n_max=0) and J2-only (n_max=2, m_max=0).
%
% Inputs:
%   r           Satellite position vector in the inertial system [m]
%   E           Transformation matrix from ECI to ECEF system
%   n_max       Maximum degree
%   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
%
% Output:
%   a           Acceleration (a=d^2r/dt^2) in the inertial system [m/s^2]
%
% Last modified:   2025/10/29   Seonwoo Son (Based on Meysam Mahooti's version)
%
%--------------------------------------------------------------------------
function a = AccelHarmonic(r, E, n_max, m_max)

global Cnm Snm C % Added global C to access C.J2 if needed

gm    = 398600.4415e9; % [m^3/s^2]; Gravitational parameter (SI)
r_ref = 6378.1363e3;   % Earth's radius [m]; Reference radius
% J2    = 0.00108263;    % J2 coefficient (use a defined value or C.J2 if available)
% If C is guaranteed to be loaded globally and contains J2:
J2 = C.J2; 

r_mag = norm(r); % Magnitude of inertial position vector

% --- Case 1: Two-Body Problem Only ---
if n_max == 2 && m_max == 0
    % Calculate Two-Body part in ECI
    a_2body = -gm / (r_mag^3) * r;

    % Calculate J2 Perturbation part using direct formula in ECEF
    r_bf = E * r; % Convert position to Body-Fixed (ECEF)
    x = r_bf(1);
    y = r_bf(2);
    z = r_bf(3);
    d = norm(r_bf); % Magnitude in ECEF (same as r_mag)
    d2 = d*d;
    d5 = d2*d2*d;
    z2_d2 = z*z / d2;

    % Direct J2 acceleration formulas in ECEF
    factor = -(3/2) * J2 * gm * r_ref^2 / d5;

    aj2x = factor * x * (1 - 5*z2_d2);
    aj2y = factor * y * (1 - 5*z2_d2);
    aj2z = factor * z * (3 - 5*z2_d2);

    a_j2_bf = [aj2x; aj2y; aj2z]; % J2 acceleration in ECEF

    % Convert J2 acceleration back to Inertial (ECI)
    a_j2_eci = E' * a_j2_bf;

    % Total acceleration is Two-Body + J2
    a = a_2body + a_j2_eci;

% --- Case 3: General Spherical Harmonics (Original Method) ---
else
    % Use the original potential gradient method for higher orders/degrees

    % Body-fixed position
    r_bf = E * r;
    % Auxiliary quantities
    d = norm(r_bf);                     % distance
    latgc = asin(r_bf(3)/d);
    lon = atan2(r_bf(2),r_bf(1));

    % Check if Legendre function needs computation
    if (n_max < 0) % Handle cases where no harmonics needed beyond 2-body (though n_max=0 covers this)
         a = -gm / (r_mag^3) * r;
         return;
    end

    [pnm, dpnm] = Legendre(n_max,m_max,latgc);
    dUdr = 0;
    dUdlatgc = 0;
    dUdlon = 0;
    q3 = 0; q2 = q3; q1 = q2;

    % Loop includes n=0, which calculates the main two-body term implicitly via C00=1
    for n=0:n_max 
        b1 = (-gm/d^2)*(r_ref/d)^n*(n+1);
        b2 =  (gm/d)*(r_ref/d)^n;
        b3 =  (gm/d)*(r_ref/d)^n;
        for m=0:m_max
            % Avoid accessing Cnm/Snm indices out of bounds if m_max < n_max
            if m > m_max 
                continue;
            end
            % Ensure Cnm and Snm have enough rows/columns up to n_max+1, m_max+1
            if size(Cnm, 1) < n+1 || size(Cnm, 2) < m+1 || size(Snm, 1) < n+1 || size(Snm, 2) < m+1
                 error('Cnm or Snm matrix size is insufficient for n_max=%d, m_max=%d', n_max, m_max);
            end

            q1 = q1 + pnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
            q2 = q2 + dpnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
            if m > 0 % Avoid calculation for m=0 where it's zero
                q3 = q3 + m*pnm(n+1,m+1)*(Snm(n+1,m+1)*cos(m*lon)-Cnm(n+1,m+1)*sin(m*lon));
            end
        end
        dUdr     = dUdr     + q1*b1;
        dUdlatgc = dUdlatgc + q2*b2;
        dUdlon   = dUdlon   + q3*b3;
        q3 = 0; q2 = q3; q1 = q2; % Reset accumulators for next n
    end

    % Body-fixed acceleration (Handle potential polar singularity)
    r2xy = r_bf(1)^2+r_bf(2)^2;
    small_threshold = 1e-16; % Threshold to avoid division by zero near poles

    if r2xy < small_threshold
        % Simplified handling near poles (may need refinement for high accuracy polar orbits)
        ax = 0; 
        ay = 0; 
        az = dUdr; % Primarily radial acceleration dominates at poles
    else
        ax = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(1)-(1/r2xy*dUdlon)*r_bf(2);
        ay = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1/r2xy*dUdlon)*r_bf(1);
        az =  1/d*dUdr*r_bf(3)+sqrt(r2xy)/d^2*dUdlatgc;
    end
    a_bf = [ax; ay; az];

    % Inertial acceleration
    a = E'*a_bf;

end % End of main if/elseif/else block

end % End of function