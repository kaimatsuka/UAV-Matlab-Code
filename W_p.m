function y = W_p(W_ENG, N_E)
% 
% DESCRIPTION: 
%   Nicolai propulsive weight estimate
%
% INPUTS:
%   W_ENG = bare engine (lbs)
%   N_E   = number of engines
%
% OUTPUTS:
%   W_p = total installed propulsion unit weight less fuel system (lb)
%
% SOURCE:
%   Nicoli weight estimate equations (MAE 154A resource)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = 2.575*(W_ENG)^0.922*N_E;

end