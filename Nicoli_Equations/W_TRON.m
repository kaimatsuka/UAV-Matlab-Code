function y = W_TRON(W_AU)
% 
% DESCRIPTION: 
%   Nicolai horizontal tail weight estimate
%
% INPUTS:
%   W_AU = bare avionics equipment weight (lbs)
%
% OUTPUTS:
%   W_TRON = total installed weight of the avionics equipment (lb)
%
% SOURCE:
%   Nicoli weight estimate equations (MAE 154A resource)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y= 2.117*(W_AU)^0.933;

end