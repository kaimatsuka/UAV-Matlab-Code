
% UAV Parameters


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
N = 6; %ultimate load factor
V_e = V_max*0.592484; %equivalent max airspeed at SL(knots)


%% UAV Parameters
%
%   These variables are hard coded UAV parameters. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_TO = 48.5; %initial weight guess of a/c (lbs)
e_H = 0.8;
C_Lmax = 1.2; %lift coefficient inital guess


% Wing --------------------------------------------------------------------

wing.b = 10; %wing span(ft)
wing.c = wing.b/13; %average wing chord length(ft)
wing.mtr = 0.12; %maximum thickness ratio
wing.t = wing.c*wing.mtr; %thickness of wing (ft)
winglet.t = 0; %thickness of winglet (ft)
wing.swp_ang = 0; %sweep angle
wing.lam = 2; %taper ratio
wing.S = wing.b*wing.c; %wing area(ft^2)(reference area)
wing.lam_q = 0; %wing quarter chord sweep
wing.x = 1.5; %(ft)
winglet.c = 0.001; %(ft)
winglet.x = 1.5; %(ft)
wing.Q = 1; %wing interference factor 
winglet.Q = 1.08;
wing.S_wet = 2.003*wing.S; %wet area for wing (ft^2)
winglet.S_wet = 0; %wet winglet area (ft^2)
A = (wing.b^2)/wing.S; %aspect ratio
K_i = pi/(A*e_H); %wing geometry

% Fuselage ----------------------------------------------------------------

fuse.L = wing.b*(3/7); %fuselage length 
fuse.W = 13/12; %fuselage max width 
fuse.D = fuse.W; %fuselage max depth
fuse.Q = 1.25; %fuselage interference factor 
fuse.r = fuse.L/fuse.W; %fuselage ratio
fuse.S_wet = pi*0.5*((fuse.L*fuse.W)+(fuse.L*fuse.W));%wet fuselage area (ft^2)

% Horizontal Tail ---------------------------------------------------------

htail.S = 0.3*wing.S; %horizontal tail area (ft)
htail.b = wing.b*0.3; %horizontal tail span (ft) 
htail.t = 0.8*wing.t; %horizontal tail max root thickness (ft)
htail.c = 0.4*wing.c; % average horizontal tail chord length (ft)
htail.x = 0.9*fuse.L; %(ft)
htail.l_T = htail.x-wing.x; %distance from wing 1/4 MAC to tail 1/4 MAC (ft)
htail.S_wet = 2.003*htail.S; %wet area for horizontal tail 
htail.lam = 0.49; %taper ratio of horizontal tail
htail.Q = 1.08; %horizontal tail interference factor 

% Vertical Tail -----------------------------------------------------------

vtail.S = htail.S/2; %vertical tail area (ft^2)
vtail.b = htail.b/2; % vertical tail span (ft)
vtail.c = 5; % average vertical tail chord length (ft)
vtail.t = 0.8*wing.t; %vertical tail max root thickness (ft)
vtail.x = 0.9*fuse.L;
vtail.Q = 1.08; %vertical tail interference factor 
vtail.S_wet = 2.003*vtail.S/2; %wet area for vertical tail
vtail.lam = 0.6; %taper ratio 
