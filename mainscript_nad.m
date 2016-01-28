%MAE 154A
%Nadia, Eugene, Anny, Kai

%NADIA: This is currently what I have for the weight estimates, I will be
%making funcQ_tions for these soon. I checked that all my values were the
%same here and on the excel spreadsheet! (01/15/2016) Functions should be
%completed by the end of the weekend, if not sooner. 

clc, clear all; close all;

% Conversion Constants

%% Load Files
load_unit_conversion
uav_params
load_enviro_parameters

%% AIRCRAFT LAYOUT

%WING----------------------------------------------------------------------

M = 0.0595760168; %Mach number

wing.Re_low = Re(rho_low,V_cruise, wing.c,mu); %Reynolds number for wing at low altitide
wing.Re_high = Re(rho_high, V_cruise,wing.c, mu); %Reynolds number for wing at high altitude
wing.Re_avg = Re(rho_avg,V_cruise,wing.c,mu); %average Reynolds number for wing
winglet.Re_avg = Re(rho_avg,V_cruise,winglet.c,mu); %average Reynolds number for winglet
wing.K = K(wing.x,wing.c,wing.t,M,wing.swp_ang); %form factor for wing
winglet.K = K(winglet.x,winglet.c,winglet.t,M,wing.swp_ang); %form factor for winglet
wing.c_f = C_f(wing.Re_avg,M); %skin friction coefficient for wing
winglet.c_f = C_f(winglet.Re_avg,M); %skin friction coefficient for winglet

%FUSELAGE------------------------------------------------------------------ 

fuse.Re_low = Re(rho_low,V_cruise,fuse.L,mu); %Reynolds number fuselage wing at low altitide
fuse.Re_high = Re(rho_high, V_cruise,fuse.L,mu); %Reynolds number for fuselage at high altitude
fuse.Re_avg = Re(rho_avg,V_cruise,fuse.L,mu); %average Reynolds number for fuselage 
fuse.K = (1+60/fuse.L^3-fuse.L/400) ; %form factor for fuselage
fuse.c_f = C_f(fuse.Re_avg,M); %skin friction coefficient for fuselage

%HORIZONTAL TAIL----------------------------------------------------------- 

htail.Re_low = Re(rho_low,V_cruise,htail.c,mu); %Reynolds number for htail at low altitide
htail.Re_high = Re(rho_high,V_cruise,htail.c, mu); %Reynolds number fot htail at high altitide
htail.Re_avg = Re(rho_avg,V_cruise,htail.c,mu); %average Reynolds number for htail
htail.K = K(htail.x,htail.c,htail.t,M,wing.swp_ang); %form factor for horizontal tail  
htail.c_f = C_f(htail.Re_avg,M); %skin friction coefficient for horizontal tail

%VERTICAL TAIL-------------------------------------------------------------

vtail.Re_low = Re(rho_low,V_cruise, vtail.c,mu);%Reynolds number fot vtail at high altitide
vtail.Re_high = Re(rho_high, V_cruise,vtail.c, mu);%Reynolds number fot vtail at low altitide
vtail.Re_avg = Re(rho_avg,V_cruise,vtail.c,mu); %average Reynolds number for vtail
vtail.K = K(vtail.x,vtail.c,vtail.t,M,wing.swp_ang); %form factor for vertical wing
vtail.c_f = C_f(vtail.Re_avg,M); %skin friction coefficient for vertical tail

%PROPULSION----------------------------------------------------------------
W_ENG = 3.5; %bare engine weight
N_E = 1; %number of engines

%ELECTRONICS/AVIONICS------------------------------------------------------ 
w_EOIR = 10; %EO/IR weight(lbs)
w_SAR = 2; %Synthetic Apperature Radar weight(lbs)
w_LiDAR = 1; %LiDAR weight(lbs)
w_ANT = 0.15; %UHF/VHF antenna weight (lbs)
w_WR = 0.25; %WaveRelay weight (lbs)
W_AU = w_EOIR+w_SAR+w_LiDAR+w_ANT+w_WR; %bare avionics equipment weight (uninstalled)

%FUEL SYSTEM---------------------------------------------------------------
F_G = 5; %total fuel in gallons
int = 23; %percent of fuel tanks that are integral 
N_T = 1; %number of separate fuel tanks

%% Drag Estimation Equations

V_drag  = V_cruise;
S_total = wing.S + htail.S;

%MISC DRAG TERMS:
C_Dmisc = 0.001;
C_DLP   = 0.001;

% Component drag parameters
K = [wing.K,fuse.K,htail.K,vtail.K,winglet.K]; %form factor vector
Q = [wing.Q,fuse.Q,htail.Q,vtail.Q,winglet.Q]; %interference factor vector
C_f = [wing.c_f,fuse.c_f,htail.c_f,vtail.c_f,winglet.c_f]; %skin friction vector
S_wet = [wing.S_wet,fuse.S_wet,htail.S_wet,vtail.S_wet,winglet.S_wet]; %wet area vector

% Compute parasite Drag
C_Dp = C_Dpi(K,Q,C_f,S_wet,wing.S,C_Dmisc,C_DLP); %parasite dragcoefficient equation

% Compute induced drag
C_l = W_TO./(0.5*rho_low*V_drag.^2*S_total);
C_Di = K_i*C_l.^2; %induced drag coefficient equation

Drag_parasite = C_Dp*0.5*rho_low*V_drag.^2*S_total;
Drag_induced  = C_Di*0.5*rho_low.*V_drag.^2*S_total;
Drag_total = Drag_parasite+Drag_induced;

%% Power Requirements

P_reqd = Drag_total*V_stall/550; %power required (HP) 

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