% calc_random_UAV.m
%
% DESCRIPTION:
%   Using the base UAV parameters, first randomly vary the adjustable 
%   parameters. Then, compute derived geometeries. 
%
% REVISION HISTORY:
%   02/10: First file is created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Wing --------------------------------------------------------------------

wing.S     = baseUAV.wing.S     + randn*wing.sd_S;     % wing area (ft^2) 
wing.A     = baseUAV.wing.A     + randn*wing.sd_A;     % aspect ratio      
wing.lam   = baseUAV.wing.lam   + randn*wing.sd_lam;   % taper ratio (must be between 0 < 1)
wing.lam_q = baseUAV.wing.lam_q + randn*wing.sd_lam_q; % wing quarter chord sweep
wing.h     = baseUAV.wing.h     + randn*wing.sd_h;     % dist from head to wing 1/4 chord (ft) TODO: where in quarter chord?
% wing.Q = 1;       % wing interference factor 
% wing.e = 0.8;     % Oswald's efficiency factor 
% %wing.gamma = 2;   % wing dihedral (deg)  GLOBAL HAWK
% 
% % airfoil properties
% %   TODO: update this once we have airfoil for wing
% wing.mtr = 0.12; % maximum thickness ratio
% wing.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
% wing.S_wet = 2.003*wing.S;  % wet area for wing (ft^2)
% 
% % Fuselage ----------------------------------------------------------------
% 
% % TODO: ass
% 
% fuse.Q = 1.25;    % fuselage interference factor TODO: cite source
% fuse.L = 4.2857;  % fuselage max length 
% fuse.W = 1.0833;  % fuselage max width 
% fuse.D = 1.0833;  % fuselage max depth
% 
% % Horizontal Tail ---------------------------------------------------------
% 
% % primary
% htail.A = 3.9;    % aspect ratio
% htail.S = 2.3077; % area (ft^2)
% htail.lam = 0.49; % taper ratio of horizontal tail (btw 0 and 1 inclusive)
% htail.lam_q = 0;  % horizontal tail sweep angle
% % htail.lam_max = 0.49; % sweep of maximum thicknes line 
% htail.h = 3.8571; % dist from head to 1/4 chord of horizontal tail (ft)
% htail.Q   = 1.08; % horizontal tail interference factor 
% 
% % airfoil properties
% %   TODO: update this once we have airfoil for horizontal tail
% htail.mtr = 0.24; % maximum thickness ratio
% htail.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
% htail.S_wet = 2.003*htail.S; % wet area for horizontal tail (ft^2)
% 
% % Vertical Tail -----------------------------------------------------------
% 
% % primay
% vtail.S = 1.1539; % area (ft^2)
% vtail.A = 1.95;   % aspect ratio (defined as b^2/S) 
% vtail.lam = 0.6;  % taper ratio 
% vtail.lam_q = 0;  % quarter chord sweep angle
% % vtail.lam_max = 0.49; % sweep of maximum thicknes line 
% vtail.Q = 1.08;   % interference factor 
% vtail.h = 3.8571; % dist from head to 1/4 chord of vertical tail (ft)
% 
% % airfoil properties
% %   TODO: update this once we have airfoil for vertical tail
% vtail.mtr = 0.96; % maximum thickness ratio
% vtail.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
% vtail.S_wet = 2.003*vtail.S;     % wet area for vertical tail (ft^2)

% Airfoil -----------------------------------------------------------------

% randomly select an airfoil
af_idx = randi(length(airfoil_options));

afoil.name  = strcat('naca',airfoil_options(af_idx));
% afoil.alpha0 = airfoils(ii).alpha0;
afoil.Cla_rad = airfoils(af_idx).Cl_alpha_rad; % (rad) 2D lift-curve slope of airfoil 
afoil.CLa_rad = afoil.Cla_rad/(1+(afoil.Cla_rad/(pi*wing.A*wing.e))); % (rad) 2D lift-curve slope of airfoil 

% % Egnine ------------------------------------------------------------------
% 
% %   no teims here
% 
% % Fuel System -------------------------------------------------------------
% 
% %   no items here
% 
% % Propeller ---------------------------------------------------------------
% 
% %   no items here
% 
% % Electronics/Payloads ----------------------------------------------------
% 
% %   no items here
% 
% % Surface Control ---------------------------------------------------------
% 
% %   no items here


%% Calculate derived geometries/variables

% Wing --------------------------------------------------------------------

wing.b = sqrt(wing.A*wing.S);     % wing span (ft)
wing.c = wing.S/wing.b;           % average wing cord (ft)   TODO: change name to wing.c_ave
wing.c_r = 2*wing.c/(1+wing.lam); % wing root chord length (ft)
wing.c_t = wing.c_r*wing.lam;     % wing tip chord length (ft)
wing.t = wing.c*wing.mtr;         % average thickness of wing (ft) TODO: change name to wing.t_ave
wing.K_i = 1/(pi*wing.A*wing.e);  % Drag coefficient

% Fuselage ----------------------------------------------------------------

fuse.r = fuse.L/fuse.W; %fuselage ratio

% arbitrarily defined 
fuse.S_wet = pi*0.5*((fuse.L*fuse.W)+(fuse.L*fuse.W));%wet fuselage area (ft^2)

% Horizontal Tail ---------------------------------------------------------

htail.b = sqrt(htail.A*htail.S); % span (ft)
htail.c = htail.S/htail.b;       % average chord length (ft) TODO: change name to htail.c_ave
htail.c_r = 2*htail.c/(1+htail.lam); % horizontal tail root chord length (ft)
htail.c_t = htail.c_r*htail.lam;     % wing tip chord length (ft)
% htail.t   = htail.c*htail.mtr;     % max thickness (average) (ft) TODO: change name to htail.t_ave
htail.t_r = htail.c_r*htail.mtr;   % max thickness at root (ft)
htail.l_T = htail.h-wing.h;      % distance from wing 1/4 MAC to tail 1/4 MAC (ft)

% Vertical Tail -----------------------------------------------------------

vtail.b = sqrt(vtail.A*vtail.S); % span (ft)
vtail.c = vtail.S/vtail.b;       % average chord length (ft)
vtail.c_r = 2*vtail.c/(1+vtail.lam); % horizontal tail root chord length (ft)
vtail.c_t = vtail.c_r*vtail.lam;     % wing tip chord length (ft)
vtail.t = vtail.c*vtail.mtr;     % max root thickness average(ft)
vtail.t_r =vtail.c_r*vtail.mtr;  % max root thickness at root (ft)

% Airfoil -----------------------------------------------------------------

% afoil.CL_max = 1.2; % Max Cl

% Egnine ------------------------------------------------------------------

% Fuel System -------------------------------------------------------------

% Propeller ---------------------------------------------------------------

% Electronics/Payloads ----------------------------------------------------

% Surface Control ---------------------------------------------------------

%% Calculate CG locations
%
% TODO: are these values correct?
%

% Adjustable CG locations
wing.x_cg  = baseUAV.wing.x_cg  + randn*wing.sd_x_cg;  % (ft)
fuse.x_cg  = baseUAV.fuse.x_cg  + randn*fuse.sd_x_cg;  % (ft)
htail.x_cg = baseUAV.htail.x_cg + randn*htail.sd_x_cg; % (ft)
vtail.x_cg = baseUAV.vtail.x_cg + randn*vtail.sd_x_cg; % (ft)
engn.x_cg  = baseUAV.engn.x_cg  + randn*engn.sd_x_cg; % (ft)
fsys.x_cg  = baseUAV.fsys.x_cg  + randn*fsys.sd_x_cg; % (ft)
prop.x_cg  = baseUAV.prop.x_cg  + randn*prop.sd_x_cg; % (ft)
payld.x_cg_EOIR  = baseUAV.payld.x_cg_EOIR  + randn*payld.sd_x_cg_EOIR;  % (ft)
payld.x_cg_SAR   = baseUAV.payld.x_cg_SAR   + randn*payld.sd_x_cg_SAR;   % (ft)
payld.x_cg_LiDAR = baseUAV.payld.x_cg_LiDAR + randn*payld.sd_x_cg_LiDAR; % (ft)
payld.x_cg_ANT   = baseUAV.payld.x_cg_ANT   + randn*payld.sd_x_cg_ANT;   % (ft)
payld.x_cg_WR    = baseUAV.payld.x_cg_WR    + randn*payld.sd_x_cg_WR;  % (ft)
payld.x_cg_IMU   = baseUAV.payld.x_cg_IMU   + randn*payld.sd_x_cg_IMU; % (ft)

% Derived CG locations
% CG location of servos
ail.x_cg_sfcl  = wing.x_cg+wing.c*(0.85-0.25)   + randn*sfcl.sd_x_cg_wing;  % (ft)
rudd.x_cg_sfcl = htail.x_cg+htail.c*(0.85-0.25) + randn*sfcl.sd_x_cg_htail; % (ft)
elev.x_cg_sfcl = vtail.x_cg+htail.c*(0.85-0.25) + randn*sfcl.sd_x_cg_vtail; % (ft)


