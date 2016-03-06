function y = W_ht(W_TO,N,l_T,S_H,b_H,t_HR)
% 
% DESCRIPTION: 
%   Nicolai horizontal tail weight estimate
%
% INPUTS:
%   W_TO  = take-off weight on wing (lbs)
%   N     = ultimate load factor (1.5 times limit load factor)
%   l_T   = distance from wing 1/4 MAC to tail 1/4 MAC (unitless)
%   S_H   = horizontal tail area (ft^2)
%   b_H   = horizontal tail span (ft)
%   t_HR  = horizontal tail max root thickness (ft) 
%
% OUTPUTS:
%   W_ht   = weight of horizontal tail (lb)
%
% SOURCE:
%   Nicoli weight estimate equations (MAE 154A resource)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ft --> in
t_HR = t_HR*12;

% NOTE: equation uses t_HR in inches
y = 127*(((((W_TO*N)/(10^5))^0.87)*((l_T/10)^0.483)*...
    ((S_H/100)^1.2)*((b_H/t_HR)^0.5))^0.458);

end
