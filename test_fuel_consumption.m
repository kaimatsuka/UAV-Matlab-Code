% test fuel consumption

clc, clear all; close all;

%% Load Files
load_unit_conversion
load_requirements
uav_params
load_enviro_parameters

%% Initial guess of take-off weight 
% W_TO = 48.5; %initial weight guess of a/c (lbs)
W_TO = 65.2; %initial weight guess of a/c (lbs)

%% Drag Calculation

% INPUT PARAMETERS for cacl_drag.m
M       = 0.0595760168; %Mach number
S_ref   = wing.S;
% rho     = rho_avg;

% calc_drag % calculate drag

%% Power Calculation

% calc_engn
engn.HP = 5.2;     %engine horse power (selected from engine)

%% Propeller Sizing

calc_propeller
P_avail = engn.HP*prop.eta_p; % Power available

%% Fuel Calculation

W_i = W_TO;
E  = [5*60 1*3600 2*3600]; %[sec]
fuel.cp = 4.2929e-07; %[1/ft] specific fuel consumption

% climb
rho     = rho_avg;
v_drag  = sqrt((76^2+132^2)/2);
calc_drag
[W_final(1) W_fuel(1)] = endu2W_fuel(W_i, prop.eta_p, E(1), fuel.cp, DRAG.C_L, DRAG.C_Dt, rho, S_ref);

% high alt scan
rho     = rho_high;
v_drag  = 132;
calc_drag
[W_final(2) W_fuel(2)] = endu2W_fuel(W_final(1), prop.eta_p, E(2), fuel.cp, DRAG.C_L, DRAG.C_Dt, rho, S_ref);


% low alt loiter
rho     = rho_low;
v_drag  = 76;
calc_drag
[W_final(3) W_fuel(3)] = endu2W_fuel(W_final(2), prop.eta_p, E(3), fuel.cp, DRAG.C_L, DRAG.C_Dt, rho, S_ref);

%sum of predicted fuel use
W_fuel_total = sum(W_fuel);