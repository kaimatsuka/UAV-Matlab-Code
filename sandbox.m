% sandbox
% nothing in this file is permanent

clear all; clc;

Wi    = 48.5;
eta_p = 0.8;
E     = 1*3600;       %[sec]
cp    = (0.7+1)/2;    %[lb/hp-hr]
cp    = cp/550/3600;  %[1/ft]
S_ref = 8.036717159;  %[ft^2] reference area

% case 1
V   = 132; %[ft/s] speed 
rho = 0.0018985; %[slug/ft^3] density at 7500 ft
CL  = 2*Wi/(rho*V^2*S_ref); %Cl at this speed
CD  = 0.04445303619; %7500 ft

W_fuel1 = endu2W_fuel(Wi, eta_p, E, cp, CL, CD, rho, S_ref);

% case 2
V   = 73.3; %[ft/s] speed
rho = 0.002309;  % density at 1000 ft
CL  = 2*Wi/(rho*V^2*S_ref); %Cl at this speed
CD  = 0.07110152435; %1000 ft

W_fuel = endu2W_fuel(Wi, eta_p, E, cp, CL, CD, rho, S_ref);