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

wing.A = 13;      % aspect ratio      
wing.S = 100/13;  % wing area (ft^2) 
wing.lam = 0.5;   % taper ratio (must be between 0 < 1)
wing.lam_q = 0;   % wing quarter chord sweep (deg)
% wing.lam_max = 0.49; % sweep of maximum thicknes line 
wing.h = 1.286;   % dist from head to wing 1/4 chord (ft) TODO: where in quarter chord?
wing.Q = 1;       % wing interference factor 
wing.e = 0.8;     % Oswald's efficiency factor 
%wing.gamma = 2;   % wing dihedral (deg)  GLOBAL HAWK

% airfoil properties
%   TODO: update this once we have airfoil for wing
wing.mtr = 0.12; % maximum thickness ratio
wing.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
wing.S_wet = 2.003*wing.S;  % wet area for wing (ft^2)

% Fuselage ----------------------------------------------------------------

fuse.Q = 1.25;    % fuselage interference factor TODO: cite source
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
htail.Q   = 1.08; % horizontal tail interference factor 

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
vtail.Q = 1.08;   % interference factor 
vtail.h = 3.8571; % dist from head to 1/4 chord of vertical tail (ft)

% airfoil properties
%   TODO: update this once we have airfoil for vertical tail
vtail.mtr = 0.12; % maximum thickness ratio
vtail.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
vtail.S_wet = 2.003*vtail.S;     % wet area for vertical tail (ft^2)

% Ailron ------------------------------------------------------------------

ail.S = 0.051*wing.S; % area (ft^2) multiplication factor ranges btwn 0 to 0.051

% Rudder ------------------------------------------------------------------

rudd.S = 0.4*wing.S; % area (ft^2) multiplication factor ranges btwn 0.3 to 0.5

% Elevator ----------------------------------------------------------------

elev.S = 0.325*wing.S; % area (ft^2) multiplication factor rangers btwn 0.3 to 0.35

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
% loadUAV.sfcl  = sfcl; % not used here

