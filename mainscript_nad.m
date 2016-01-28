%MAE 154A
%Nadia, Eugene, Anny, Kai

%NADIA: This is currently what I have for the weight estimates, I will be
%making funcQ_tions for these soon. I checked that all my values were the
%same here and on the excel spreadsheet! (01/15/2016) Functions should be
%completed by the end of the weekend, if not sooner. 

clc, clear all; close all;


%% Load Files
load_unit_conversion
uav_params
load_enviro_parameters

%% Drag Calculation

% INPUT PARAMETERS for cacl_drag.m
M       = 0.0595760168;   %Mach number
v_drag  = V_cruise;
S_total = wing.S + htail.S;
rho     = rho_avg;

calc_drag % calculate drag

%% WEIGHT CALCULATION------------------------------------------------------

%PROPULSION----------------------------------------------------------------
W_ENG = 3.5; %bare engine weight
N_E = 1; %number of engines

%ELECTRONICS/AVIONICS------------------------------------------------------

w_EOIR  = 10; %EO/IR weight(lbs)
w_SAR   = 2; %Synthetic Apperature Radar weight(lbs)
w_LiDAR = 1; %LiDAR weight(lbs)
w_ANT   = 0.15; %UHF/VHF antenna weight (lbs)
w_WR    = 0.25; %WaveRelay weight (lbs)
W_AU    = w_EOIR+w_SAR+w_LiDAR+w_ANT+w_WR; %bare avionics equipment weight (uninstalled)

%FUEL SYSTEM---------------------------------------------------------------
F_G = 5; %total fuel in gallons
int = 23; %percent of fuel tanks that are integral 
N_T = 1; %number of separate fuel tanks

%% Nicolai Estimate Weight Equations

% wing = W_w(lam, S_W, W_TO, N, A, lam_q, mtr, V_e); %wing 
% 
% fuselage = W_f(W_TO,N,L,W,D,V_e); %fuselage
% 
% horizontal_tail = W_ht(W_TO,N,l_T,S_H,b_H,t_HR); %horizontal tail
% 
% vertical_tail = W_vt(W_TO,N,S_V,b_V,t_VR); %vertical tail
% 
% propulsion = W_p(W_ENG,N_E); %propulsion
% 
% avionics = W_TRON(W_AU); %electronics/avionics
% 
% fuel_system = W_FS(F_G,int,N_T, N_E); %fuel system
% 
% surface_controls = 1.08*(W_TO^0.7); %surface controls(powered)
% 
% electrical_system = W_ES(fuel_system,avionics); %electrical system

% T = table(weight, fuselage, horizontal_tail, vertical_tail, propulsion,...
%     avionics, fuel_system, surface_controls, electrical_system, 'RowNames',Value)