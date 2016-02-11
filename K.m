function y = K(mtl,mtr,M,swp_max)
% 
% DESCRIPTION:
%   This calculates form factor for subsonic drag estimation for wing,
%   tail,strut, and pylon.
%
% INPUT:
%   mtl     = chordwise location of airfoil max thickness location
%   mtr     = maximum thickness ratio
%   M       = mach number 
%   swp_max = sweep angle of max thickness line
%
% REFERENCE:
%   Raymer, pg 435 (chap 12.5.4)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = (1+0.6/mtl*mtr+100*mtr^4)*...
    (1.34*M.^0.18*cos(swp_max)^0.28);

end