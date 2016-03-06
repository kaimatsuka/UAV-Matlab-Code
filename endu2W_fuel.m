function [W_final W_fuel] = endu2W_fuel(Wi, eta_p, E, cp, CL, CD, rho, S_ref)
% endu2mfuel.m
%
% Description:
%   Calculate the fuel weight from the endurance and initial weight.
%
% INPUTS:
%   Wi    = initial weight(lb)
%   eta_p = propeller efficiency
%   E     = endurance (sec)
%   cp    = specific heat of fuel
%   CL    = lift coefficient
%   CD    = drag coefficient
%   rho   = air density (slug/ft^3)
%   S_ref = reference area (ft^2)
%
% OUTPUTS:
%   W_fuel = weight of fuel (lb)
%

W_final = (Wi^-0.5+E*cp*CD/eta_p/CL^1.5/sqrt(2*rho*S_ref))^-2;
W_fuel = Wi-W_final;

end

