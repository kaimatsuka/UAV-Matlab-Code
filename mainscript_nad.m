%MAE 154A
%Nadia, Eugene, Anny, Kai
%
%NADIA: This is currently what I have for the weight estimates, I will be
%making funcQ_tions for these soon. I checked that all my values were the
%same here and on the excel spreadsheet! (01/15/2016) Functions should be
%completed by the end of the weekend, if not sooner. 
%
% REVISION HISTORY:
%   01/29: fixed unit converstion errors in Nicoli (htail, vtail
%          thickness)
%

clc, clear all; close all;

%% Load Files
load_unit_conversion
load_requirements
uav_params
load_enviro_parameters

%% Initial guess of take-off weight 
W_TO = 48.5; %initial weight guess of a/c (lbs)

%% Drag Calculation

% INPUT PARAMETERS for cacl_drag.m
M       = 0.0595760168;   %Mach number
v_drag  = V_cruise;
S_total = wing.S + htail.S;
rho     = rho_avg;

calc_drag % calculate drag

%% Power Calculation

% calc_engn
engn.HP = 5.2;     %engine horse power (selected from engine)

%% Propeller Sizing

calc_propeller
P_avail = engn.HP*prop.eta_p; % Power available

%% Weight Calculation

% Engine(s) ---------------------------------------------------------------

engn.W_bare = 3.5; % bare engine weight
engn.N_E = 1;      % number of engines

%ELECTRONICS/AVIONICS------------------------------------------------------

w_EOIR  = 10; %EO/IR weight(lbs)
w_SAR   = 2;  %Synthetic Apperature Radar weight(lbs)
w_LiDAR = 1;  %LiDAR weight(lbs)
w_ANT   = 0.15; %UHF/VHF antenna weight (lbs)
w_WR    = 0.25; %WaveRelay weight (lbs)
W_AU    = w_EOIR+w_SAR+w_LiDAR+w_ANT+w_WR; % bare avionics equipment weight (uninstalled)

%FUEL SYSTEM---------------------------------------------------------------

F_G = 0.668403; %total fuel (ft^3) TODO: calculate this from endurance
int = 23; %percent of fuel tanks that are integral TODO: find out how much is integral
N_T = 1; %number of separate fuel tanks 

%% Nicolai Estimate Weight Equations

WEIGHT.wing  = W_w(wing.lam, wing.S, W_TO, N, wing.A, wing.lam_q, wing.mtr, V_e); %wing 
WEIGHT.fuse  = W_f(W_TO,N,fuse.L,fuse.W,fuse.D,V_e); %fuselage
WEIGHT.htail = W_ht(W_TO,N, htail.l_T, htail.S, htail.b, htail.t); %horizontal tail
WEIGHT.vtail = W_vt(W_TO,N,vtail.S,vtail.b,vtail.t); %vertical tail
WEIGHT.engn  = W_p(engn.W_bare,engn.N_E); %propulsion
WEIGHT.avion = W_TRON(W_AU); %electronics/avionics
WEIGHT.fsys  = W_FS(F_G,int,N_T, engn.N_E); %fuel system
WEIGHT.sc    = W_SC(W_TO); %surface controls(powered)
% WEIGHT.esys  = W_ES(WEIGHT.fsys,WEIGHT.avion); %electrical system

%TODO: WEIGHT.fuel

weight_vec = [WEIGHT.wing  WEIGHT.fuse   WEIGHT.htail  WEIGHT.vtail ...
              WEIGHT.engn  WEIGHT.avion  WEIGHT.fsys   WEIGHT.sc];
weight_labels = {'wing',   'fuse',     'htail',       'vtail', ...
                'engine', 'avionics', 'fuel system', 'surface controls'};
WEIGHT.total = sum(weight_vec);

figure()
pie(weight_vec, weight_labels); hold on;

%TODO:
% WEIGHT.fuel = fuel carried by UAV 
% WEIGHT.gear = gears needed for launching and catching (optional)

