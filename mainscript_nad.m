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
% uav_params
% load_UAV_parameters
load_base_UAV
load_variation_parameters
calc_random_UAV
load_enviro_parameters

%% Initial guess of take-off weight 
% W_TO = 48.5; %initial weight guess of a/c (lbs)
W_TO = 33.9; %initial weight guess of a/c (lbs)

%% Drag Calculation

% INPUT PARAMETERS for cacl_drag.m
M       = 0.0595760168;   %Mach number
v_drag  = V_cruise;
S_ref   = wing.S;
rho     = rho_avg;

calc_drag % calculate drag

%% Weight Calculation
calc_weight_estimate

weight_vec = [WEIGHT.wing  WEIGHT.fuse   WEIGHT.htail  WEIGHT.vtail ...
              WEIGHT.engn  WEIGHT.avion  WEIGHT.fsys   WEIGHT.sc ...
              WEIGHT.fuel  WEIGHT.prop];
weight_labels = {'wing',  'fuselage',         'horizontal tail', 'vertical tail', ...
                'engine', 'payload/avionics', 'fuel system',     'surface controls',...
                'fuel',   'propeller'};

WEIGHT.total2 = sum(weight_vec);

if 1
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
            10 ...
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