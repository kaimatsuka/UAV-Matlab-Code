% load_base_UAV.m
%
% DESCRIPTION:
%   This file generate and define the basic parameters of base UAV. 
%
% DEPENDENCY:
%   Before calling this file, call
%       load_airfoils
%       engine_directory
%
% REVISION HISTORY:
%   02/10: First file is created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Wing --------------------------------------------------------------------

wing.S = 15;      % wing area (ft^2) 
wing.A = 10;      % aspect ratio      
wing.lam = 0.5;   % taper ratio (must be between 0 < 1)
wing.lam_q = 0;   % wing quarter chord sweep (deg)
% wing.lam_max = 0.49; % sweep of maximum thicknes line 
wing.h_q = 1.286+1.5;   % dist from head to wing 1/4 chord at root (ft)

%wing.gamma = 2;   % wing dihedral (deg)  GLOBAL HAWK

% airfoil properties
airfoilw.ind = 1; % base airfoil is 1!

% Fuselage ----------------------------------------------------------------

fuse.L = 4.6315+3;  % fuselage max length 
fuse.W = 1.0833;  % fuselage max width 
fuse.D = 1.0833;  % fuselage max depth

% Horizontal Tail ---------------------------------------------------------

% primary
htail.A = 3.9;    % aspect ratio
htail.S = 2.3077+5; % area (ft^2)
htail.lam = 1; % taper ratio of horizontal tail (btw 0 and 1 inclusive)
htail.lam_q = 0;  % horizontal tail sweep angle (deg)
% htail.lam_max = 0.49; % sweep of maximum thicknes line 
htail.h = 3.8571+2.5; % dist from head to 1/4 chord of horizontal tail (ft)

aifoilh.ind = 1; % default is 1!

% Vertical Tail -----------------------------------------------------------

% primay
vtail.S = 1.1539; % area (ft^2)
vtail.A = 1.95;   % aspect ratio (defined as b^2/S) 
vtail.lam = 1;  % taper ratio 
vtail.lam_q = 0;  % quarter chord sweep angle (deg)
% vtail.lam_max = 0.49; % sweep of maximum thicknes line 
vtail.h = 3.8571+3; % dist from head to 1/4 chord of vertical tail (ft)

airfoilv.ind = 1;

% Fuel --------------------------------------------------------------------

fuel.W = 6; %[lb] fuel weight (calculated using test_fuel.m file)

% Fuel System -------------------------------------------------------------

%   no items here

% Engine ------------------------------------------------------------------

engn.ind = 1; % choose 1 (out of 9) as base engine

% Propeller ---------------------------------------------------------------

prop.D = 1.62;   % propeller diameter (ft)
prop.pitch = 15; % propeller pitch (ft)

% Electronics/Payloads ----------------------------------------------------

% LiDAR, ANT, WR, and IMU are movable payload, while EOIR and SAR are not.
payld.x_cg_LiDAR = 0.075*fuse.L; % (ft)
payld.x_cg_ANT   = 0.05*fuse.L;  % (ft)
payld.x_cg_WR    = 0.05*fuse.L;  % (ft)
payld.x_cg_IMU   = 0.05*fuse.L; % (ft)

% Aileron -----------------------------------------------------------------

sfcl.ail.S = 1.1250; % area (ft^2) 3%~12% of wing area
sfcl.ail.c = 0.2449; % aileron chord (ft) 15%~30% of wing chord
sfcl.ail.maxallow = 30; % degrees
sfcl.ail.minallow = -25; % degrees

% Rudder ------------------------------------------------------------------

sfcl.rudd.S = 0.2596; % area (ft^2) 15%~35% of vertical tail area
sfcl.rudd.c = 0.2115; % rudder chord (ft) 15%~40% of vertical tail chord
sfcl.rudd.maxallow = 24; % degrees
sfcl.rudd.minallow = -20; % degrees

% Elevator ----------------------------------------------------------------

sfcl.elev.S = 0.6346; % area (ft^2) 15%~40% of horizontal tail area
sfcl.elev.c = 0.2308; % elevator chord (ft) 20%~40% of horizontal tail chord
sfcl.elev.maxallow = 25; % degrees
sfcl.elev.minallow = -20; % degrees

%% calculate derived parameters
calc_non_tunable_parameters

%% Create baseUAV structure ------------------------------------------------

baseUAV.wing  = wing;
baseUAV.fuse  = fuse;
baseUAV.htail = htail;
baseUAV.vtail = vtail;
baseUAV.fuse  = fuse;
baseUAV.engn  = engn;
baseUAV.fsys  = fsys;
baseUAV.fuel  = fuel;
baseUAV.prop  = prop;
baseUAV.payld = payld;
baseUAV.sfcl  = sfcl;

% clear wing fuse htail vtail fuse engn fsys prop payld