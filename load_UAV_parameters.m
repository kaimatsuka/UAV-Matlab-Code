% load_UAV_params.m
%
% DESCRIPTION:
%    Hard coded UAV parameters and their derived parameters
%
% REVISION HISTORY:
%   01/29: moved W_TO to main_script, fix typo for K_i equation. 
%   02/09: chordwise location of the airofil max thickness loactions are 
%          defined for wing, vtail, htail
%   02/10: separated adjustable variables from derived variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Standard deviation/variation range of variation
%
% These are assigned to parameters that designer can vary.
%
% 
% wing.sd_A     = 0;
% wing.sd_S     = 0;
% wing.sd_lam   = 0;
% wing.sd_lam_q = 0;
% wing.sd_mtr   = 0;
% wing.x        = 0;

%% Adjustable Variables

% Wing --------------------------------------------------------------------

wing.A = 13;      % aspect ratio      
wing.S = 100/13;  % wing area (ft^2) 
wing.lam = 0.5;   % taper ratio (must be between 0 < 1)
wing.lam_q = 0;   % wing quarter chord sweep
% wing.lam_max = 0.49; % sweep of maximum thicknes line 
wing.x = 1.286;   % dist from head to wing 1/4 chord (ft) TODO: where in quarter chord?
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
fuse.L = 4.2857;  % fuselage max length 
fuse.W = 1.0833;  % fuselage max width 
fuse.D = 1.0833;  % fuselage max depth

% Horizontal Tail ---------------------------------------------------------

% primary
htail.A = 3.9;    % aspect ratio
htail.S = 2.3077; % area (ft^2)
htail.lam = 0.49; % taper ratio of horizontal tail (btw 0 and 1 inclusive)
htail.lam_q = 0;  % horizontal tail sweep angle
% htail.lam_max = 0.49; % sweep of maximum thicknes line 
htail.x = 3.8571; % dist from head to 1/4 chord of horizontal tail (ft)
htail.Q   = 1.08; % horizontal tail interference factor 

% airfoil properties
%   TODO: update this once we have airfoil for horizontal tail
htail.mtr = 0.24; % maximum thickness ratio
htail.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
htail.S_wet = 2.003*htail.S; % wet area for horizontal tail (ft^2)

% Vertical Tail -----------------------------------------------------------

% primay
vtail.S = 1.1539; % area (ft^2)
vtail.A = 1.95;   % aspect ratio (defined as b^2/S) 
vtail.lam = 0.6;  % taper ratio 
vtail.lam_q = 0;  % quarter chord sweep angle
% vtail.lam_max = 0.49; % sweep of maximum thicknes line 
vtail.Q = 1.08;   % interference factor 
vtail.x = 3.8571; % dist from head to 1/4 chord of vertical tail (ft)

% airfoil properties
%   TODO: update this once we have airfoil for vertical tail
vtail.mtr = 0.96; % maximum thickness ratio
vtail.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
vtail.S_wet = 2.003*vtail.S;     % wet area for vertical tail (ft^2)

% Airfoil -----------------------------------------------------------------

%   no teims here

% Egnine ------------------------------------------------------------------

%   no items here

% Fuel System -------------------------------------------------------------

%   no items here

% Propeller ---------------------------------------------------------------

%   no items here

% Electronics/Payloads ----------------------------------------------------

%   no items here

% Surface Control ---------------------------------------------------------

%   no items here

%% Derived geometries/variables

% Wing --------------------------------------------------------------------

wing.b = sqrt(wing.A*wing.S);     % wing span (ft)
wing.c = wing.S/wing.b;           % average wing cord (ft)   TODO: change name to wing.c_ave
wing.c_r = 2*wing.c/(1+wing.lam); % wing root chord length (ft)
wing.c_t = wing.c_r*wing.lam;     % wing tip chord length (ft)
wing.t = wing.c*wing.mtr;         % average thickness of wing (ft) TODO: change name to wing.t_ave
wing.K_i = 1/(pi*wing.A*wing.e);  % Drag coefficient
wing.x_cg = wing.x;               % Wing CG location wrt nose (ft) TODO: is this correct?

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
htail.l_T = htail.x-wing.x;      % distance from wing 1/4 MAC to tail 1/4 MAC (ft)
htail.x_cg = htail.x;            % tail CG location (ft) TODO: is this correct?

% Vertical Tail -----------------------------------------------------------

vtail.b = sqrt(vtail.A*vtail.S); % span (ft)
vtail.c = vtail.S/vtail.b;       % average chord length (ft)
vtail.c_r = 2*vtail.c/(1+vtail.lam); % horizontal tail root chord length (ft)
vtail.c_t = vtail.c_r*vtail.lam;     % wing tip chord length (ft)
vtail.t = vtail.c*vtail.mtr;     % max root thickness average(ft)
vtail.t_r =vtail.c_r*vtail.mtr;  % max root thickness at root (ft)
vtail.x_cg = vtail.x;            % vertical tail CG location (ft) TODO: is this correct?

% Airfoil -----------------------------------------------------------------

% afoil.CL_max = 1.2; % Max Cl

% Egnine ------------------------------------------------------------------

engn.x_cg = 0.95*fuse.L; % engine CG location (assume 95% of fuselage)

% Fuel System -------------------------------------------------------------

fsys.x_cg = 0.65*fuse.L; % fuel system CG location (ft)

% Propeller ---------------------------------------------------------------

prop.x_cg = (1+0.025)*fuse.L; % propeller CG location (ft)

% Electronics/Payloads ----------------------------------------------------

% CG locations
payld.x_cg_EOIR  = 0.025*fuse.L; % EO/IR CG location (ft)
payld.x_cg_SAR   = 0.05*fuse.L;  % Synthetic Apperature Radar CG location (ft)
payld.x_cg_LiDAR = 0.075*fuse.L; % LiDAR CG location (ft)
payld.x_cg_ANT   = 0.05*fuse.L;  % UHF/VHF antenna location (ft)
payld.x_cg_WR    = 0.05*fuse.L;  % WaveRelay CG location (ft)
payld.x_cg_IMU   = 0.05*fuse.L;  % IMU CG location (ft)

% Surface Control ---------------------------------------------------------

sfcl.x_cg_wing  = wing.x_cg+wing.c*(0.85-0.25); % assume 85% of average wing chord
sfcl.x_cg_htail = htail.x_cg+htail.c*(0.85-0.25); % assume 85% of average wing chord
sfcl.x_cg_vtail = vtail.x_cg+htail.c*(0.85-0.25); % assume 85% of average wing chord



