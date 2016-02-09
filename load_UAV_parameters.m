% load_UAV_params.m
%
% DESCRIPTION:
%    Hard coded UAV parameters and their derived parameters
%
% REVISION HISTORY:
%   01/29: moved W_TO to main_script, fix typo for K_i equation. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Standard deviation/variation range of variation
%
% These are assigned to parameters that designer can vary.
%

wing.sd_A     = 0;
wing.sd_S     = 0;
wing.sd_lam   = 0;
wing.sd_lam_q = 0;
wing.sd_mtr   = 0;
wing.x        = 0;

%% Calculate aircraft geometries

% Wing --------------------------------------------------------------------

% primary
wing.A = 13;      %aspect ratio         SCAN EAGLE
wing.S = 100/13;  %wing area (ft^2)     SCAN EAGLE
wing.lam = 0.5;   %taper ratio (must be between 0 < 1)
wing.lam_q = 0;   %wing quarter chord sweep
wing.mtr = 0.12;  %maximum thickness ratio
wing.x = 1.286;   %dist from head to wing 1/4 chord (ft)
wing.Q = 1;       %wing interference factor 
wing.e = 0.8;     %Oswald's efficiency factor 
%wing.gamma = 2;   %wing dihedral (deg)  GLOBAL HAWK

% derived
wing.b = sqrt(wing.A*wing.S);     % wing span (ft)
wing.c = wing.S/wing.b;           % average wing cord (ft)
wing.c_r = 2*wing.c/(1+wing.lam); % wing root chord length (ft)
wing.c_t = wing.c_r*wing.lam;     % wing tip chord length (ft)
wing.t = wing.c*wing.mtr;         % thickness of wing (ft)
wing.S_wet = 2.003*wing.S;        % wet area for wing (ft^2)
wing.K_i = 1/(pi*wing.A*wing.e);  % Drag coefficient
wing.x_cg = wing.x;               % Wing CG location wrt nose (ft)

% Winglet -----------------------------------------------------------------
% 
% % primary
% winglet.Q = 1.08;   %winglet interference factor TODO: where to get this?
% winglet.c = 0.001;  %winglet cord length (ft)
% winglet.t = 0;      %thickness of winglet (ft)
% winglet.S_wet = 0;  %wet winglet area (ft^2) TODO: put some value
% winglet.lam_q = 0;  %sweep angle of winglet
% 
% % derived
% winglet.x = wing.x + 0; %dist from head to winglet 1/4 chord (ft)

% Fuselage ----------------------------------------------------------------

% primary
fuse.Q = 1.25;  %fuselage interference factor TODO: cite source
fuse.L = 4.2857;  %fuselage length 
fuse.W = 1.0833;  %fuselage max width 

% derived
fuse.D = fuse.W;  %fuselage diameter TODO: is this same thing as fuse.W?
fuse.r = fuse.L/fuse.W; %fuselage ratio

% arbitrarily selected
fuse.S_wet = pi*0.5*((fuse.L*fuse.W)+(fuse.L*fuse.W));%wet fuselage area (ft^2)
fuse.x_cg = fuse.L/2;   %fuselage CG location wrt nose (ft)

% Horizontal Tail ---------------------------------------------------------

% primary
htail.lam = 0.49; %taper ratio of horizontal tail (btw 0 and 1 inclusive)
htail.Q   = 1.08; %horizontal tail interference factor 
htail.lam_q = 0;  %horizontal tail sweep angle
htail.S = 2.3077; % area (ft)

% derived
htail.b = wing.b*0.3; % span (ft) 
htail.t = 0.8*wing.t; % max root thickness (ft)
htail.c = 0.4*wing.c; % average chord length (ft)
htail.x = 0.9*fuse.L; % dist from head to 1/4 chord of horizontal tail (ft)
htail.l_T = htail.x-wing.x; % distance from wing 1/4 MAC to tail 1/4 MAC (ft)
htail.S_wet = 2.003*htail.S;% wet area for horizontal tail 
htail.x_cg = htail.x; % horizontal tail CG location (ft)

% arbitrariliy selected


% Vertical Tail -----------------------------------------------------------

% primay
vtail.c = 5;     % average chord length (ft)
vtail.lam = 0.6; % taper ratio 
vtail.Q = 1.08;  % interference factor 
vtail.lam_q = 0; % quarter chord sweep angle

% derived
vtail.S = htail.S/2;    % area (ft^2)
vtail.b = htail.b/2;    % span (ft)
vtail.t = 0.8*wing.t;   % max root thickness (ft)
vtail.x = 0.9*fuse.L;   % dist from head to 1/4 chord of vertical tail (ft)
vtail.S_wet = 2.003*vtail.S/2; % wet area for vertical tail
vtail.x_cg = vtail.x;  %  vertical tail CG location (ft)

% Airfoil -----------------------------------------------------------------

afoil.CL_max = 1.2; % Max Cl

% Egnine ------------------------------------------------------------------

engn.W_bare = 3.5; % bare engine weight
engn.N_E = 1;      % number of engines
engn.x_cg = 0.95*fuse.L; % engine CG location (assume 95% of fuselage)

% Fuel System -------------------------------------------------------------

fsys.x_cg = 0.65*fuse.L; % fuel system CG location (ft)

% Propeller ---------------------------------------------------------------

prop.x_cg = (1+0.025)*fuse.L; % propeller CG location (ft)

% Electronics/Payloads ----------------------------------------------------

% weights
payld.w_EOIR  = 12.57; % EO/IR weight(lbs)
payld.w_SAR   = 2;     % Synthetic Apperature Radar weight(lbs)
payld.w_LiDAR = 1;     % LiDAR weight(lbs)
payld.w_ANT   = 0.15;  % UHF/VHF antenna weight (lbs)
payld.w_WR    = 0.25;  % WaveRelay weight (lbs)
payld.w_IMU   = 0.11;  % Ellipse E mini (lbs)

% CG locations
payld.x_cg_EOIR  = 0.025*fuse.L; % EO/IR CG location (ft)
payld.x_cg_SAR   = 0.05*fuse.L;  % Synthetic Apperature Radar CG location (ft)
payld.x_cg_LiDAR = 0.075*fuse.L; % LiDAR CG location (ft)
payld.x_cg_ANT   = 0.05*fuse.L;  % UHF/VHF antenna location (ft)
payld.x_cg_WR    = 0.05*fuse.L;  % WaveRelay CG location (ft)
payld.x_cg_IMU   = 0.05*fuse.L;  % IMU CG location (ft)

% Pyaload/Avionics total weight
payld.w_total = payld.w_EOIR + payld.w_SAR + payld.w_LiDAR + payld.w_ANT ...
                + payld.w_WR + payld.w_IMU;

% Surface Control ---------------------------------------------------------

sfcl.x_cg_wing = wing.x_cg+wing.c*(0.85-0.25); % assume 85% of average wing chord
sfcl.x_cg_htail = htail.x_cg+htail.c*(0.85-0.25); % assume 85% of average wing chord
sfcl.x_cg_vtail = vtail.x_cg+htail.c*(0.85-0.25); % assume 85% of average wing chord