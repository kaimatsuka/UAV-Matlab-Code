% test_drag.m
%
%   Stand alone test code that will plot drag and power required as a 
%   fucnction of velocity.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all;

load_UAV_parameters;
load_unit_conversion;

%% Input Parameter
V = [50:10:200];      % specify velocity range for plotting
S_total = S_W + S_H;

%WING----------------------------------------------------------------------

M      = 0.0595760168; %Mach number

Re_Wlow = Re(rho_low,V_cruise, c_W,mu); %Reynolds number for wing at low altitide
Re_Whigh = Re(rho_high, V_cruise,c_W, mu); %Reynolds number for wing at high altitude
Re_Wavg = Re(rho_avg,V_cruise,c_W,mu); %average Reynolds number for wing
Re_WWavg = Re(rho_avg,V_cruise,c_winglet,mu); %average Reynolds number for winglet
K_wing = K(x_wing,c_W,t_W,M,swp_ang); %form factor for wing
K_winglet = K(x_winglet,c_winglet,t_winglet,M,swp_ang); %form factor for winglet
C_fW = C_f(Re_Wavg,M); %skin friction coefficient for wing
C_fWW = C_f(Re_WWavg,M); %skin friction coefficient for winglet

%FUSELAGE------------------------------------------------------------------ 

Re_Flow = Re(rho_low,V_cruise, L,mu); %Reynolds number fuselage wing at low altitide
Re_Fhigh = Re(rho_high, V_cruise,L, mu); %Reynolds number for fuselage at high altitude
Re_Favg = Re(rho_avg,V_cruise,L,mu); %average Reynolds number for fuselage 
K_fuse = (1+60/L^3-L/400) ; %form factor for fuselage
C_fF = C_f(Re_Favg,M); %skin friction coefficient for fuselage

%HORIZONTAL TAIL----------------------------------------------------------- 

Re_Hlow = Re(rho_low,V_cruise, c_H,mu); %Reynolds number for htail at low altitide
Re_Hhigh = Re(rho_high, V_cruise,c_H, mu); %Reynolds number fot htail at high altitide
Re_Havg = Re(rho_avg,V_cruise,c_H,mu); %average Reynolds number for htail
K_H = K(x_Htail,c_H,t_HR,M,swp_ang); %form factor for horizontal tail  
C_fH = C_f(Re_Havg,M); %skin friction coefficient for horizontal tail

%VERTICAL TAIL-------------------------------------------------------------

Re_Vlow = Re(rho_low,V_cruise, c_V,mu);%Reynolds number fot vtail at high altitide
Re_Vhigh = Re(rho_high, V_cruise,c_V, mu);%Reynolds number fot vtail at low altitide
Re_Vavg = Re(rho_avg,V_cruise,c_H,mu); %average Reynolds number for vtail
K_V = K(x_Vtail,c_H,t_HR,M,swp_ang); %form factor for vertical wing
C_fV = C_f(Re_Vavg,M); %skin friction coefficient for vertical tail

%PROPULSION----------------------------------------------------------------
W_ENG = 3.5; %bare engine weight
N_E = 1; %number of engines



%% Drag Estimation Equations

%MISC DRAG TERMS:
C_Dmisc = 0.001;
C_DLP   = 0.001;

% Component drag parameters
K     = [K_wing,K_fuse,K_H,K_V,K_winglet]; %form factor vector
Q     = [Q_wing, Q_fuse,Q_H,Q_V,Q_winglet]; %interference factor vector
C_f   = [C_fW,C_fF,C_fH,C_fV,C_fWW]; %skin friction vector
S_wet = [S_Wwet,S_Fwet,S_Hwet,S_Vwet,S_WWwet]; %wet area vector

% Compute parasite Drag
C_Dp = C_Dpi(K,Q,C_f,S_wet,S_W,C_Dmisc,C_DLP); %parasite dragcoefficient equation

% Compute induced drag
C_l = W_A./(0.5*rho_low*V.^2*S_total);
C_Di = K_i*C_l.^2; %induced drag coefficient equation

Drag_parasite = C_Dp*0.5*rho_low*V.^2*S_total;
Drag_induced  = C_Di*0.5*rho_low.*V.^2*S_total;
Drag_total = Drag_parasite+Drag_induced;

%% Power Required
Power_req_parasite = Drag_parasite.*V*lbfts2hp;
Power_req_induced  = Drag_induced.*V*lbfts2hp;
Power_req_total    = Drag_total.*V*lbfts2hp;

%% Plot Drag
figure()
plot(V,Drag_parasite, 'b'); hold on; grid on;
plot(V,Drag_induced, 'r');
plot(V,Drag_total, 'color', [0 0.5 0]);
xlabel('Velocity(ft/s)'),ylabel('Drag(lb)');
legend('Parasite','Induced','Total Drag');
title('Drag vs. Velcoity');

figure()
plot(V,Power_req_parasite, 'b'); hold on; grid on;
plot(V,Power_req_induced, 'r'); 
plot(V,Power_req_total, 'color', [0 0.5 0]); 
xlabel('Velocity(ft/s)'),ylabel('Power(HP)');
legend('Parasite','Induced','Total Power Req');
title('Power Required vs. Velcoity');