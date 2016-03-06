% load_requiremetns.m
%
% DESCRIPTION:
%    UAV design requirement parameters that are derived from mission
%    objectives.
%
% REVISION HISTORY:
%   01/29: This file is created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_unit_conversion

% primary requirements
V_stall  = 45*mph2fps; % stall speed(ft/s)
V_loiter = 50*mph2fps; % loiter speed(ft/s)
V_cruise = 90*mph2fps; % cruise speed(ft/s)
RC_max   = 15;         % maximum rate of climb (ft/s)
V_max    = 95*mph2fps; % max speed(ft/s)
Ran      = 100;    % range (ft) TODO: fix this number
E_max    = 1*3600; % endurance at 7500 ft(seconds) TODO: fix this number
E_min    = 2*3600; % endurance at 1000 ft(seconds) TODO: fix this number
LF_V = 8; % load factor at V
N = 6;    % ultimate load factor (limit load is 4)


% derived requirements
E_T = E_max+E_min; % total endurance(seconds)
V_e = V_max;       % equivalent max airspeed at SL(ft/s)
V_launch = V_stall*1.15; % Launch speed (ft/s)
N_lim = N/1.5; % limit load factor