function y = W_f(W_TO,N,L,W,D,V_e)
% 
% DESCRIPTION: 
%   Nicolai fuselage weight estimate
%
% INPUTS:
%   W_TO  = take-off weight on wing (lbs)
%   N     = ultimate load factor (1.5 times limit load factor)
%   L     = fuselage length (ft)
%   W     = fuselage max width (ft)
%   D     = fuselage max depth (ft)
%   V_e   = equivalent max airspeed at SL (ft/s)
%
% OUTPUTS:
%   W_f   = weight of fuselage (lb)
%
% SOURCE:
%   Nicoli weight estimate equations (MAE 154A resource)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ft/s --> knots
V_e = V_e*0.592484;

% NOTES: equation below use V_e in knots
y = 200*(((((W_TO*N)/(10^5))^0.286)*((L/10)^0.857)*...
    ((W+D)/10)*((V_e/100)^0.338))^1.1);

end