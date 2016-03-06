function y = W_w(lam, S_W, W_TO, N, A, lam_q, mtr, V_e)
% 
% DESCRIPTION: 
%   Nicolai wing weight estimate
%
% INPUTS:
%   lam   = wing taper ratio
%   S_W   = wing area (ft^2)
%   W_TO  = take-off weight on wing (lbs)
%   N     = ultimate load factor (1.5 times limit load factor)
%   A     = wing aspect ratio
%   lam_q = wing quarter chord sweep angle (degree)
%   mtr   = maximum thickness ratio
%   V_e   = equivalent max airspeed at SL (ft/s)
%
% OUTPUTS:
%   W_w   = weight of wing (lb)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ft/s --> knots
V_e = V_e*0.592484;

% NOTE: V_e below should be in knots 
y = 96.948*(((((W_TO*N)/(10^5))^0.65)*((A/cosd(lam_q))^0.57)*...
((S_W/100)^0.61)*(((1+lam)/(2*mtr))^0.36)*((1+(V_e/500))^0.5))^0.993);

end