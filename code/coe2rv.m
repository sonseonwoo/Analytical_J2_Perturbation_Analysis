% ------------------------------------------------------------------------------
%
%                           function coe2rv
%
%  this function finds the position and velocity vectors in geocentric
%    equatorial (ijk) system given the classical orbit elements.
%
%  author        : david vallado                  
%  revisions
%    vallado     - add constant file use          29 jun 2003
%    modified    - use global struct C            09 sep 2025
%
%  inputs          description                    range / units
%    p           - semilatus rectum               km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to pi rad
%    omega       - longitude of ascending node    0.0  to 2pi rad
%    argp        - argument of perigee            0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
%    truelon     - true longitude            (ce) 0.0  to 2pi rad
%    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
%
%  outputs       : [6 x 1]
%    r           - ijk position vector            km
%    v           - ijk velocity vector            km / s
%
% ------------------------------------------------------------------------------
function [r,v] = coe2rv ( p,ecc,incl,omega,argp,nu,arglat,truelon,lonper )

    % ---- use global constants ----
    global C
    mu = C.GM_Earth * 1e-9;   % [m^3/s^2] → [km^3/s^2]

    small = 1e-10;            % tolerance for zero checks

    % -------------------------------------------------------------
    %       determine what type of orbit is involved and set up
    %       angles for the special cases.
    % -------------------------------------------------------------
    if ( ecc < small )
        % ----------------  circular equatorial  ------------------
        if (incl<small) || ( abs(incl-pi)< small )
            argp = 0.0;
            omega= 0.0;
            nu   = truelon;
        else
            % --------------  circular inclined  ------------------
            argp= 0.0;
            nu  = arglat;
        end
    else
        % ---------------  elliptical equatorial  -----------------
        if ( ( incl<small) || (abs(incl-pi)<small) )
            argp = lonper;
            omega= 0.0;
        end
    end

    % ----------  form pqw position and velocity vectors ----------
    cosnu= cos(nu);
    sinnu= sin(nu);
    temp = p / (1.0  + ecc*cosnu);

    rpqw = [temp*cosnu;
            temp*sinnu;
            0.0];

    if ( abs(p) < 1.0e-4 )
        p= 1.0e-4;
    end

    vpqw = [ -sinnu*sqrt(mu/p);
              (ecc + cosnu)*sqrt(mu/p);
               0.0 ];

    % ----------------  perform transformation to ijk  ------------
    

    R1 = @(theta) [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
    R2 = @(theta) [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
    R3 = @(theta) [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
    tempvec = R3(-argp)*rpqw;
    tempvec = R1(-incl)*tempvec;
    r       = R3(-omega)*tempvec;
    
    tempvec = R3(-argp)*vpqw;
    tempvec = R1(-incl)*tempvec;
    v       = R3(-omega)*tempvec;

end
