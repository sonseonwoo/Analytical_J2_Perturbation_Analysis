%  - - - - - - -
%   i a u A n p
%  - - - - - - -
%
%  Normalize angle into the range 0 <= a < 2pi.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     a        double     angle (radians)
%
%  Returned (function value):
%              double     angle in range 0-2pi
%
%  This revision:  2021 May 11
%
%  SOFA release 2021-05-12
%
%  Copyright (C) 2021 IAU SOFA Board.  See notes at end.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = iauAnp(a)

global C

w = rem(a, C.D2PI);
if (w < 0)
    w = w + C.D2PI;
end

