%MAE 154A
%Nadia, Eugene, Anny, Kai

%NADIA: This is currently what I have for the weight estimates, I will be
%making funcQ_tions for these soon. I checked that all my values were the
%same here and on the excel spreadsheet! (01/15/2016) Functions should be
%completed by the end of the weekend, if not sooner. 

clc, clear all;
%% PARAMETERS

e = 0.8;
W_A = 48.5; %initial weight guess of a/c (lbs)
C_Lmax = 1.2; %lift coefficient inital guess
V_stall = 45*5280/3600; %stall speed(ft/s)
V_loiter = 50*5280/3600; %loiter speed(ft/s)
V_cruise = 90*5280/3600; %cruise speed(ft/s)
V_max = 95*5280/3600; %max speed(ft/s)
Ran = 100; %range (ft)
E_max = 1*3600; %endurance at 7500 ft(seconds)
E_min = 2*3600; %endurance at 1000 ft(seconds)
E_T = E_max+E_min; %total endurance(seconds)
LF_V = 8; %load factor at V
rho_low = 0.002309; %low altitude density (slugs/ft^3)
rho_high = 0.0018985; %low altitude density (slugs/ft^3)
rho_avg = (rho_low+rho_high)/2; %average density (slugs/ft^3)
gamma = 1.4;
temp = 50; %temperature (fahrenheit)
R = 1716.59; %gas constant (ft^2/s^2*R)
a_s = 1107.815608; %speed of sound (ft/s)
M = 0.0595760168; %Mach number
mu = 3.82*10^-7; %dynamic viscosity of air 

%% AIRCRAFT LAYOUT

%WING----------------------------------------------------------------------
b_W = 10; %wing span(ft)
c_W = b_W/13; %average wing chord length(ft)
mtr = 0.12; %maximum thickness ratio
t_W = c_W*mtr; %thickness of wing (ft)
t_WW = 0; %thickness of winglet (ft)
swp_ang = 0; %sweep angle
lam = 2; %taper ratio
S_W = (2*W_A)/(rho_low*C_Lmax*V_stall^2); %wing area(ft^2)(reference area)
W_TO = 48.5; %take off weight
N = 6; %ultimate load factor
A = (b_W^2)/S_W; %aspect ratio
lam_q = 0; %wing quarter chord
V_e = V_max*0.592484; %equivalent max airspeed at SL(knots) 
x_wing = 1.5; %(ft)
c_winglet = 0.001; %(ft)
x_winglet = 1.5; %(ft)
t_winglet = 0;
Q_wing = 1; %wing interference factor 
Q_winglet = 1.08;
S_Wwet = 2.003*S_W; %wet area for wing (ft^2)
S_WWwet = 0; %wet winglet area (ft^2)
Re_Wlow = Re(rho_low,V_cruise, c_W,mu); %Reynolds number for wing at low altitide
Re_Whigh = Re(rho_high, V_cruise,c_W, mu); %Reynolds number for wing at high altitude
Re_Wavg = Re(rho_avg,V_cruise,c_W,mu); %average Reynolds number for wing
Re_WWavg = Re(rho_avg,V_cruise,c_winglet,mu); %average Reynolds number for winglet
K_wing = K(x_wing,c_W,t_W,M,swp_ang); %form factor for wing
K_winglet = K(x_winglet,c_winglet,t_winglet,M,swp_ang); %form factor for winglet
C_fW = C_f(Re_Wavg,M); %skin friction coefficient for wing
C_fWW = C_f(Re_WWavg,M); %skin friction coefficient for winglet

%FUSELAGE------------------------------------------------------------------ 
L = b_W*(3/7); %fuselage length 
W = 13/12; %fuselage max width 
D = W; %fuselage max depth
Q_fuse = 1.25; %fuselage interference factor 
f_r = L/W; %fuselage ratio
S_Fwet = pi*0.5*((L*W)+(L*W));%wet fuselage area (ft^2)
Re_Flow = Re(rho_low,V_cruise, L,mu); %Reynolds number fuselage wing at low altitide
Re_Fhigh = Re(rho_high, V_cruise,L, mu); %Reynolds number for fuselage at high altitude
Re_Favg = Re(rho_avg,V_cruise,L,mu); %average Reynolds number for fuselage 
K_fuse = (1+60/L^3-L/400) ; %form factor for fuselage
C_fF = C_f(Re_Favg,M); %skin friction coefficient for fuselage

%HORIZONTAL TAIL----------------------------------------------------------- 
S_H = 0.3*S_W; %horizontal tail area (ft)
b_H = b_W*0.3; %horizontal tail span (ft) 
t_HR = 0.8*t_W; %horizontal tail max root thickness (ft)
c_H = 0.4*c_W; % average horizontal tail chord length (ft)
x_Htail = 0.9*L; %(ft)
l_T = x_Htail-x_wing; %distance from wing 1/4 MAC to tail 1/4 MAC (ft)
K_i = pi/(A*e); %wing geometry
Q_H = 1.08; %horizontal tail interference factor 
S_Hwet = 2.003*S_H; %wet area for horizontal tail 
lam_H = 0.49; %taper ratio of horizontal tail
Re_Hlow = Re(rho_low,V_cruise, c_H,mu); %Reynolds number for htail at low altitide
Re_Hhigh = Re(rho_high, V_cruise,c_H, mu); %Reynolds number fot htail at high altitide
Re_Havg = Re(rho_avg,V_cruise,c_H,mu); %average Reynolds number for htail
K_H = K(x_Htail,c_H,t_HR,M,swp_ang); %form factor for horizontal tail  
C_fH = C_f(Re_Havg,M); %skin friction coefficient for horizontal tail

%VERTICAL TAIL-------------------------------------------------------------
S_V = S_H/2; %vertical tail area (ft^2)
b_V = b_H/2; % vertical tail span (ft)
c_V = 5; % average vertical tail chord length (ft)
t_VR = 0.8*t_W; %vertical tail max root thickness (ft)
x_Vtail = 0.9*L;
Q_V = 1.08; %vertical tail interference factor 
S_Vwet = 2.003*S_V/2; %wet area for vertical tail
lam_V = 0.6; %taper ratio 
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

Drag = C_D*0.5*rho_low*V_stall^2*S_W; %Total drag (lbf)

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