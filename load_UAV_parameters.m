% UAV Parameters

%% Environment Parameter
%
%   These parameters are environmental parameters, independent of aircraft
%   or mission requirements.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gamma = 1.4; %specific heat ratio
rho_low  = 0.002309; %low altitude density (slugs/ft^3) 1000 ft
rho_high = 0.0018985; %low altitude density (slugs/ft^3) 7500 ft
rho_avg  = (rho_low+rho_high)/2; %average density (slugs/ft^3)
R        = 1716.59; %gas constant (ft^2/s^2*R)
a_s      = 1107.815608; %speed of sound (ft/s) 50F
temp     = 50; %temperature (fahrenheit)
mu       = 3.82*10^-7; %dynamic viscosity of air 


%% Mission Derived Parameters
%   These variables are derived directly from mission requirements

V_stall = 45*5280/3600; %stall speed(ft/s)
V_loiter = 50*5280/3600; %loiter speed(ft/s)
V_cruise = 90*5280/3600; %cruise speed(ft/s)
V_max = 95*5280/3600; %max speed(ft/s)
Ran = 100; %range (ft)
E_max = 1*3600; %endurance at 7500 ft(seconds)
E_min = 2*3600; %endurance at 1000 ft(seconds)
E_T = E_max+E_min; %total endurance(seconds)
LF_V = 8; %load factor at V



%% UAV Parameters
%
%   These variables are hard coded UAV parameters. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_A    = 48.5; %initial weight guess of a/c (lbs)
e_H    = 0.8;
C_Lmax = 1.2; %lift coefficient inital guess


% Wing --------------------------------------------------------------------

%   Temporal note on how to use matlab structure
%         WING.b = 10;              % wing span(ft)
%         WING.c = WING.b/13;       % average wing chord length(ft)
%         WING.mtr = 0.12;          % maximum thickness ratio
%         WING.t = WING.c*WING.mtr; % thickness of wing (ft)
%         WING.t_WW = 0;            % thickness of winglet (ft)
%         ....

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


% Fuselage ----------------------------------------------------------------

L = b_W*(3/7); %fuselage length 
W = 13/12; %fuselage max width 
D = W; %fuselage max depth
Q_fuse = 1.25; %fuselage interference factor 
f_r = L/W; %fuselage ratio
S_Fwet = pi*0.5*((L*W)+(L*W));%wet fuselage area (ft^2)

% Horizontal Tail ---------------------------------------------------------

S_H = 0.3*S_W; %horizontal tail area (ft)
b_H = b_W*0.3; %horizontal tail span (ft) 
t_HR = 0.8*t_W; %horizontal tail max root thickness (ft)
c_H = 0.4*c_W; % average horizontal tail chord length (ft)
x_Htail = 0.9*L; %(ft)
l_T = x_Htail-x_wing; %distance from wing 1/4 MAC to tail 1/4 MAC (ft)
K_i = pi/(A*e_H); %wing geometry
S_Hwet = 2.003*S_H; %wet area for horizontal tail 
lam_H = 0.49; %taper ratio of horizontal tail
Q_H = 1.08; %horizontal tail interference factor 

% Vertical Tail -----------------------------------------------------------

S_V = S_H/2; %vertical tail area (ft^2)
b_V = b_H/2; % vertical tail span (ft)
c_V = 5; % average vertical tail chord length (ft)
t_VR = 0.8*t_W; %vertical tail max root thickness (ft)
x_Vtail = 0.9*L;
Q_V = 1.08; %vertical tail interference factor 
S_Vwet = 2.003*S_V/2; %wet area for vertical tail
lam_V = 0.6; %taper ratio 
