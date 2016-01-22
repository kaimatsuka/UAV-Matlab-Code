%MAE 154A
%Nadia, Eugene, Anny, Kai

%NADIA: This is currently what I have for the weight estimates, I will be
%making funcQ_tions for these soon. I checked that all my values were the
%same here and on the excel spreadsheet! (01/15/2016) Functions should be
%completed by the end of the weekend, if not sooner. 

clc, clear all;

%% Load Files
UAV_parameters

%% AIRCRAFT LAYOUT

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

%MISC DRAG TERMS:
C_Dmisc = 0.001;
C_DLP = 0.001;

%--------------------------------------------------------------------------
%PARASITE DRAG EQUATIONs:
 
K = [K_wing,K_fuse,K_H,K_V,K_winglet]; %form factor vector

Q = [Q_wing, Q_fuse,Q_H,Q_V,Q_winglet]; %interference factor vector

C_f = [C_fW,C_fF,C_fH,C_fV,C_fWW]; %skin friction vector
 
S_wet = [S_Wwet,S_Fwet,S_Hwet,S_Vwet,S_WWwet]; %wet area vector

C_Dp = C_Dpi(K,Q,C_f,S_wet,S_W,C_Dmisc,C_DLP); %parasite dragcoefficient equation

C_Di = K_i*C_Lmax^2; %induced drag coefficient equation

C_D = C_Dp+C_Di; %Total drag coefficient 

Drag = C_D*0.5*rho_low*V_stall^2*S_W; %Total drag at low altitude, (lbf)

%% Power Requirements

P_reqd = Drag*V_stall/550; %power required (HP) 

%% Nicolai Estimate Weight Equations

wing = W_w(lam, S_W, W_TO, N, A, lam_q, mtr, V_e); %wing 

fuselage = W_f(W_TO,N,L,W,D,V_e); %fuselage

horizontal_tail = W_ht(W_TO,N,l_T,S_H,b_H,t_HR); %horizontal tail

vertical_tail = W_vt(W_TO,N,S_V,b_V,t_VR); %vertical tail

propulsion = W_p(W_ENG,N_E); %propulsion

avionics = W_TRON(W_AU); %electronics/avionics

fuel_system = W_FS(F_G,int,N_T, N_E); %fuel system

surface_controls = 1.08*(W_TO^0.7); %surface controls(powered)

electrical_system = W_ES(fuel_system,avionics); %electrical system

% T = table(weight, fuselage, horizontal_tail, vertical_tail, propulsion,...
%     avionics, fuel_system, surface_controls, electrical_system, 'RowNames',Value)