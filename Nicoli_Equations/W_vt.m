function y = W_vt(W_TO,N,S_V,b_V,t_VR)
% 
% DESCRIPTION: 
%   Nicolai vertical tail weight estimate
%
% INPUTS:
%   W_TO  = take-off weight on wing (lbs)
%   N     = ultimate load factor (1.5 times limit load factor)
%   S_V   = vertical tail area (ft^2)
%   b_V   = vertical tail span (ft)
%   t_VR  = vertical tail max root thickness (ft) 
%
% OUTPUTS:
%   W_vt   = weight of vertical tail (lb)
%
% SOURCE:
%   Nicoli weight estimate equations (MAE 154A resource)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ft --> inches
t_VR = t_VR*12;

% NOTE: equation below uses t_VR in inches
y = 98.5*(((((W_TO*N)/(10^5))^0.87)*...
    ((S_V/100)^1.2)*((b_V/t_VR)^0.5))^0.458);

end