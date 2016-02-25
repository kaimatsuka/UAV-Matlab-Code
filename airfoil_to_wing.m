function [airfoils] = airfoil_to_wing(airfoils, AR, ep )
% 
% DESCRIPTION: 
%   Converts 2D airfoil coefficients and statistics into 3D wing
%   characteristics
%
% INPUTS:
%   AR           = aspect ratio of the wing
%   ep           = wing efficiency factor 
%
% OUTPUTS: (within airfoils structure)
%   CL_alpha_deg = lift-curve slope of the wing (in degrees)
%   CL_alpha_rad = life-curve slope of the wing (in radians)
%   CL           = CL profile for all alphas
%   CLmax        = calculates max CL value
%   CL_maxClCd   = determines CL of minimum drag condition
%   alpha_maxClCd= determines alpha of minimum drag condition
%
% SOURCE:
%   MAE 154S lift-curve slope calculations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:length(airfoils)
    % Calculate the slope both /deg and /rad --> Cl_alpha
    airfoils(ii).CL_alpha_rad = airfoils(ii).Cl_alpha_rad/(1+(airfoils(ii).Cl_alpha_rad/(pi*AR*ep)));
    airfoils(ii).CL_alpha_deg = airfoils(ii).CL_alpha_rad*(pi/180);
    
    %Calculate the 3D CL values
    airfoils(ii).CL = airfoils(ii).CL_alpha_deg*(airfoils(ii).alpha-airfoils(ii).alpha0);
    
    % Calculate the max CL value (stall CL)
    airfoils(ii).CLmax = max(airfoils(ii).CL);

    % Find CL and alpha corresponding to max Cl/Cd value
    idx = find(airfoils(ii).ClCd - airfoils(ii).maxClCd == 0);
    airfoils(ii).CL_maxClCd = airfoils(ii).CL(idx);
    airfoils(ii).alpha_maxClCd = airfoils(ii).alpha(idx);

end

