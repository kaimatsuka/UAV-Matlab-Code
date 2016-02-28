% load_base_UAV.m
%
% DESCRIPTION:
%   This file generate and define the basic parameters of base UAV. 
%
% REVISION HISTORY:
%   02/10: First file is created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Wing --------------------------------------------------------------------

wing.A = 10;      % aspect ratio      
wing.S = 15;  % wing area (ft^2) 
wing.lam = 0.5;   % taper ratio (must be between 0 < 1)
wing.lam_q = 0;   % wing quarter chord sweep (deg)
% wing.lam_max = 0.49; % sweep of maximum thicknes line 
wing.h = 1.286;   % dist from head to wing 1/4 chord (ft) TODO: where in quarter chord?
wing.e = 0.8;     % Oswald's efficiency factor 
%wing.gamma = 2;   % wing dihedral (deg)  GLOBAL HAWK

% airfoil properties
%   TODO: update this once we have airfoil for wing
wing.mtr = 0.12; % maximum thickness ratio
wing.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
wing.S_wet = 2.003*wing.S;  % wet area for wing (ft^2)

% Fuselage ----------------------------------------------------------------

% fuse.L = 4.2857;  % fuselage max length 
fuse.L = 4.6315;  % fuselage max length 
fuse.W = 1.0833;  % fuselage max width 
fuse.D = 1.0833;  % fuselage max depth

% Horizontal Tail ---------------------------------------------------------

% primary
htail.A = 3.9;    % aspect ratio
htail.S = 2.3077; % area (ft^2)
htail.lam = 0.49; % taper ratio of horizontal tail (btw 0 and 1 inclusive)
htail.lam_q = 0;  % horizontal tail sweep angle (deg)
% htail.lam_max = 0.49; % sweep of maximum thicknes line 
htail.h = 3.8571; % dist from head to 1/4 chord of horizontal tail (ft)

% airfoil properties
%   TODO: update this once we have airfoil for horizontal tail
htail.mtr = 0.12; % maximum thickness ratio
htail.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
htail.S_wet = 2.003*htail.S; % wet area for horizontal tail (ft^2)

% Vertical Tail -----------------------------------------------------------

% primay
vtail.S = 1.1539; % area (ft^2)
vtail.A = 1.95;   % aspect ratio (defined as b^2/S) 
vtail.lam = 0.6;  % taper ratio 
vtail.lam_q = 0;  % quarter chord sweep angle (deg)
% vtail.lam_max = 0.49; % sweep of maximum thicknes line 
vtail.h = 3.8571; % dist from head to 1/4 chord of vertical tail (ft)

% airfoil properties
%   TODO: update this once we have airfoil for vertical tail
vtail.mtr = 0.12; % maximum thickness ratio
vtail.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
vtail.S_wet = 2.003*vtail.S;     % wet area for vertical tail (ft^2)

% Propeller ---------------------------------------------------------------

prop.h = 4.6315; % propeller location
prop.D = 1.62; % propeller diameter (ft)

% Ailron ------------------------------------------------------------------

ail.S = 0.051*wing.S; % area (ft^2) multiplication factor ranges btwn 0 to 0.051

% Rudder ------------------------------------------------------------------

rudd.S = 0.4*wing.S; % area (ft^2) multiplication factor ranges btwn 0.3 to 0.5

% Elevator ----------------------------------------------------------------

elev.S = 0.325*wing.S; % area (ft^2) multiplication factor rangers btwn 0.3 to 0.35

%% Calculate derived geometries/variables

% Wing --------------------------------------------------------------------

wing.b = sqrt(wing.A*wing.S);     % wing span (ft)
wing.c = wing.S/wing.b;           % average wing cord (ft)   TODO: change name to wing.c_ave
wing.c_r = 2*wing.c/(1+wing.lam); % wing root chord length (ft)
wing.c_t = wing.c_r*wing.lam;     % wing tip chord length (ft)
wing.t_r = wing.c_r*wing.mtr;     % max thickness at root (ft)
wing.t_t = wing.c_t*wing.mtr;     % max thickness at tip (ft)
wing.t = wing.c*wing.mtr;         % average thickness of wing (ft) TODO: change name to wing.t_ave
wing.K_i = 1/(pi*wing.A*wing.e);  % Drag coefficient
wing.x_cg = wing.h;               % Wing CG location wrt nose (ft) TODO: is this correct?

wing.h_q_t = wing.h+tand(wing.lam_q)*wing.b/2;  % wing c/4 tip location (ft)
wing.h_l_r = wing.h-wing.c_r/4;                 % wing leading edge root location (ft)
wing.h_l_t = wing.h_q_t-wing.c_t/4;             % wing leading edge tip location (ft)
wing.lam_l = atand((wing.h_l_t-wing.h_l_r)/(wing.b/2)); % wing leading edge sweep location (deg)

% Fuselage ----------------------------------------------------------------

fuse.r = fuse.L/fuse.W; %fuselage ratio

% arbitrarily defined 
fuse.S_wet = pi*0.5*((fuse.L*fuse.W)+(fuse.L*fuse.W));%wet fuselage area (ft^2)
fuse.x_cg  = fuse.L/2;   %CG location wrt nose (ft) TODO: is this correct?

% Horizontal Tail ---------------------------------------------------------

htail.b = sqrt(htail.A*htail.S); % span (ft)
htail.c = htail.S/htail.b;       % average chord length (ft) TODO: change name to htail.c_ave
htail.c_r = 2*htail.c/(1+htail.lam); % horizontal tail root chord length (ft)
htail.c_t = htail.c_r*htail.lam;     % wing tip chord length (ft)
% htail.t   = htail.c*htail.mtr;     % max thickness (average) (ft) TODO: change name to htail.t_ave
htail.t_r = htail.c_r*htail.mtr;   % max thickness at root (ft)
htail.t_t = htail.c_t*htail.mtr;   % max thickness ratio (ft)s
htail.l_T = htail.h-wing.h;      % distance from wing 1/4 MAC to tail 1/4 MAC (ft)
htail.x_cg = htail.h;            % tail CG location (ft) TODO: is this correct?

htail.h_q_t = htail.h+tand(htail.lam_q)*htail.b/2;
htail.h_l_r = htail.h    -htail.c_r/4;
htail.h_l_t = htail.h_q_t-htail.c_t/4;
htail.lam_l = atand((htail.h_l_t-htail.h_l_r)/(htail.b/2));

% Vertical Tail -----------------------------------------------------------

vtail.b = sqrt(vtail.A*vtail.S); % span (ft)
vtail.c = vtail.S/vtail.b;       % average chord length (ft)
vtail.c_r = 2*vtail.c/(1+vtail.lam); % horizontal tail root chord length (ft)
vtail.c_t = vtail.c_r*vtail.lam;     % wing tip chord length (ft)
vtail.t = vtail.c*vtail.mtr;     % max root thickness average(ft)
vtail.t_r =vtail.c_r*vtail.mtr;  % max root thickness at root (ft)
vtail.t_t = vtail.c_t*vtail.mtr; % max root thickness at tip (ft)
vtail.x_cg = vtail.h;            % vertical tail CG location (ft) TODO: is this correct?

vtail.h_q_t = vtail.h+tand(vtail.lam_q)*vtail.b;
vtail.h_l_r = vtail.h    -vtail.c_r/4;
vtail.h_l_t = vtail.h_q_t-vtail.c_t/4;
vtail.lam_l = atand((vtail.h_l_t-vtail.h_l_r)/(vtail.b/2));

% Airfoil -----------------------------------------------------------------

% afoil.CL_max = 1.2; % Max Cl

% Engine ------------------------------------------------------------------

engn.x_cg = 0.95*fuse.L; % engine CG location (assume 95% of fuselage)

% Fuel System -------------------------------------------------------------

fsys.x_cg = 0.65*fuse.L; % fuel system CG location (ft)

% Propeller ---------------------------------------------------------------

prop.h    = fuse.L; % beginning of prop is located at end of fuselge
prop.x_cg = (1+0.025)*fuse.L; % propeller CG location (ft)

% CG locations ------------------------------------------------------------

wing.x_cg  = wing.h;   % (ft)
fuse.x_cg  = fuse.L/2; % (ft)
htail.x_cg = htail.h;  % (ft)
vtail.x_cg = vtail.h;  % (ft)
engn.x_cg = 0.95*fuse.L;  % (ft)
fsys.x_cg = 0.65*fuse.L;  % (ft)
prop.x_cg = 1.025*fuse.L; % (ft)
payld.x_cg_EOIR  = 0.025*fuse.L; % (ft)
payld.x_cg_SAR   = 0.05*fuse.L;  % (ft)
payld.x_cg_LiDAR = 0.075*fuse.L; % (ft)
payld.x_cg_ANT   = 0.05*fuse.L;  % (ft)
payld.x_cg_WR    = 0.05*fuse.L;  % (ft)
payld.x_cg_IMU   = 0.05*fuse.L; % (ft)

% Surface control CG is dervied from geometries of wing surface and 
% control surfaces. Calculated in calc_random_UAV
% sfcl.x_cg_wing  = wing.x_cg+wing.c*(0.85-0.25); % assume 85% of average wing chord
% sfcl.x_cg_htail = htail.x_cg+htail.c*(0.85-0.25); % assume 85% of average htail chord
% sfcl.x_cg_vtail = vtail.x_cg+htail.c*(0.85-0.25); % assume 85% of average htail chord

% Create baseUAV structure ------------------------------------------------

baseUAV.wing  = wing;
baseUAV.fuse  = fuse;
baseUAV.htail = htail;
baseUAV.vtail = vtail;
baseUAV.fuse  = fuse;
baseUAV.engn  = engn;
baseUAV.fsys  = fsys;
baseUAV.prop  = prop;
baseUAV.payld = payld;
% baseUAV.sfcl  = sfcl; % not used here

% clear wing fuse htail vtail fuse engn fsys prop payld