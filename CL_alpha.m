function [CL_alpha_wing_rad,CL_alpha_wing_deg ] = CL_alpha( Cl_alpha_rad, AR, ep )
% 
% DESCRIPTION: 
%   Determine CL_alpha based off of Cl_alpha of airfoil, aspect ratio, and
%   efficiency of wing
%
% INPUTS:
%   Cl_alpha_rad = lift-curve slope of the airfoil (in radians)
%   AR           = aspect ratio of the wing
%   ep           = wing efficiency factor 
%
% OUTPUTS:
%   CL_alpha_deg = lift-curve slope of the wing (in degrees)
%   CL_alpha_rad = life-curve slope of the wing (in radians)
%
% SOURCE:
%   MAE 154S lift-curve slope calculations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CL_alpha_wing_rad = Cl_alpha_rad/(1+(Cl_alpha_rad/(pi*AR*ep)));
CL_alpha_wing_deg = CL_alpha_wing_rad*(pi/180);
end

