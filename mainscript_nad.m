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
% W_TO = 48.5; %initial weight guess of a/c (lbs)
W_TO = 39.74; %initial weight guess of a/c (lbs)

%% Drag Calculation

% INPUT PARAMETERS for cacl_drag.m
M       = 0.0595760168;   %Mach number
v_drag  = V_cruise;
S_ref   = wing.S;
rho     = rho_avg;

calc_drag % calculate drag

%% Power Calculation

% calc_engn
engn.HP = 5.2;     %engine horse power (selected from engine)

%% Propeller Sizing

calc_propeller
P_avail = engn.HP*prop.eta_p; % Power available

%% Fuel Calculation

fuel.cp = 0.85/550/3600;  %[1/ft] specific fuel consumption
fuel.rho = 6.073; %[lbm/gallon] density of fuel for octane gas
fuel.W = 3.5;             %[lb] fuel weight (calculated using test_fuel.m file)
fuel.V = fuel.W/fuel.rho; %[gallon] volume of fuel
fuel.V = fuel.V*gallon2ft3; %[ft^3] volume of fuel

%% Fuel Calculation

fsys.int = 0.00; %percent of fuel tanks that are integral TODO: find out how much is integral
fsys.N_T = 1;    %number of separate fuel tanks

%% Weight Estimation

% Nicoli Weight Estimat Equation
WEIGHT.wing  = W_w(wing.lam, wing.S, W_TO, N, wing.A, wing.lam_q, wing.mtr, V_e); %wing 
WEIGHT.fuse  = W_f(W_TO,N,fuse.L,fuse.W,fuse.D,V_e); %fuselage
WEIGHT.htail = W_ht(W_TO,N,htail.l_T, htail.S, htail.b, htail.t); %horizontal tail
WEIGHT.vtail = W_vt(W_TO,N,vtail.S,vtail.b,vtail.t); %vertical tail
WEIGHT.fsys  = W_FS(fuel.V,fsys.int,fsys.N_T, engn.N_E); %fuel system %TODO: update with actual fuel tank 
% WEIGHT.engn  = W_p(engn.W_bare,engn.N_E); %propulsion
% WEIGHT.avion = W_TRON(W_AU); %electronics/avionics
% WEIGHT.sc    = W_SC(W_TO); %surface controls (powered)
% WEIGHT.esys  = W_ES(WEIGHT.fsys,WEIGHT.avion); %electrical system

% Sum of components 
WEIGHT.engn  = 1.1*engn.W_bare*engn.N_E; % 10% mounting fudge factor
WEIGHT.avion = payld.w_total;  % electronics/avionics
WEIGHT.sc    = 0.4915*W_TO^(2/3); % Nadia Equation %TODO: update with actual servos
WEIGHT.fuel  = fuel.W;
WEIGHT.prop  = 0.75; %TODO: incorporate equation


weight_vec = [WEIGHT.wing  WEIGHT.fuse   WEIGHT.htail  WEIGHT.vtail ...
              WEIGHT.engn  WEIGHT.avion  WEIGHT.fsys   WEIGHT.sc ...
              WEIGHT.fuel  WEIGHT.prop];
weight_labels = {'wing',  'fuselage',         'horizontal tail', 'vertical tail', ...
                'engine', 'payload/avionics', 'fuel system',     'surface controls',...
                'fuel',   'propeller'};

% Detailed weight list
w_detail_vec = [WEIGHT.wing ...
                WEIGHT.fsys ...
                WEIGHT.htail ...
                WEIGHT.vtail ...
                WEIGHT.engn ...
                WEIGHT.fsys ...
                WEIGHT.prop ... 
                payld.w_EOIR ...
                payld.w_SAR ...
                payld.w_LiDAR ...
                payld.w_ANT ...
                payld.w_IMU ...
                payld.w_WR ...
                WEIGHT.sc*0.75 ...
                WEIGHT.sc*0.125 ...
                WEIGHT.sc*0.125];

WEIGHT.total = sum(weight_vec);

if 0
    figure()
    pie(weight_vec, weight_labels); hold on;
end

%% Volume Calculation

% TODO: 

%% CG Calculation

% total_cg = sum(component_weight*component_cg)/sum(component_ewight)

% TODO
x_cg_vec = [wing.x_cg ...
            fuse.x_cg ... 
            htail.x_cg ... 
            vtail.x_cg ... 
            engn.x_cg ... 
            fsys.x_cg ... 
            prop.x_cg ...
            payld.x_cg_EOIR ...
            payld.x_cg_SAR ...
            payld.x_cg_LiDAR ...
            payld.x_cg_ANT ...
            payld.x_cg_IMU ...
            payld.x_cg_WR ...
            sfcl.x_cg_wing ...
            sfcl.x_cg_htail ...
            sfcl.x_cg_vtail];
        
        
total_cg = x_cg_vec*w_detail_vec'/sum(w_detail_vec);