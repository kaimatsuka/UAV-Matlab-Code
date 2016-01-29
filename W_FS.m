function y = W_FS(F_vol,int,N_T, N_E) 
% 
% DESCRIPTION: 
%   Nicolai fuel system weight estimate. This includes pumps, lines and
%   tanks.
%
% INPUTS:
%   F_vol = max fuel volume (ft^3) 
%   int   = percent of fuel tanks that are integral (eg 0.2 if 20%)
%   N_T   = number of separate fuel tanks
%   N_E   = number of engines
%
% OUTPUTS:
%   W_FS   = weight of horizontal tail (lb)
%
% SOURCE:
%   Nicoli weight estimate equations (MAE 154A resource)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ft^3 --> gallon
F_vol = F_vol*7.48052;

% NOTE: equation below uses F_vol in gallon
y = 2.49*(((F_vol^0.6)*((1/(1+int))^0.3)*(N_T^0.2)*(N_E^0.13))^1.21);

end