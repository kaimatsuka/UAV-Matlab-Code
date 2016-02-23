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

wing.S     = baseUAV.wing.S     + (rand-0.5)*wing.sd_S; % wing area (ft^2) 
wing.A     = baseUAV.wing.A     + (rand-0.5)*wing.sd_A; % aspect ratio      
wing.lam   = baseUAV.wing.lam   + (rand-0.5)*wing.sd_lam; % taper ratio (must be between 0 < 1)
wing.lam_q = baseUAV.wing.lam_q + (rand-0.5)*wing.sd_lam_q; % wing quarter chord sweep
wing.h     = baseUAV.wing.h     + (rand-0.5)*wing.sd_h; % dist from head to wing 1/4 chord (ft) TODO: where in quarter chord?
wing.Q     = baseUAV.wing.Q     + (rand-0.5)*wing.sd_Q; % wing interference factor 
wing.e = 0.8;     % Oswald's efficiency factor 
%wing.gamma = 2;   % wing dihedral (deg)  GLOBAL HAWK

% airfoil properties
airf_index_w       = randi(12); % randomly select airfoil
airfoilw.CL_alpha = airfoils(airf_index_w).CL_alpha_deg; % CL_alpha (/deg)
airfoilw.CLmax    = airfoils(airf_index_w).CLmax; % CLmax (wing)
airfoilw.maxthick = airfoils(airf_index_w).GEO.max_thick; % max thickness (relative to chord length)
airfoilw.maxthick_loc = airfoils(airf_index_w).GEO.max_thick_location; % max thickness location (relative to chord length)
airfoilw.perim    = airfoils(airf_index_w).GEO.perimeter; % circumference/perimeter of airfoil

wing.mtr = 0.12; % maximum thickness ratio
wing.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
wing.S_wet = 2.003*wing.S;  % wet area for wing (ft^2)

% Fuselage ----------------------------------------------------------------

fuse.Q = baseUAV.fuse.Q  + (rand-0.5)*fuse.sd_Q;    % fuselage interference factor TODO: cite source
fuse.L = baseUAV.fuse.L  + (rand-0.5)*fuse.sd_L;  % fuselage max length 
fuse.W = baseUAV.fuse.W  + (rand-0.5)*fuse.sd_W;  % fuselage max width 
fuse.D = baseUAV.fuse.D  + (rand-0.5)*fuse.sd_D;  % fuselage max depth

% Horizontal Tail ---------------------------------------------------------

% primary
htail.A = baseUAV.htail.A              + (rand-0.5)*htail.sd_A;    % aspect ratio
htail.S = baseUAV.htail.S              + (rand-0.5)*htail.sd_S; % area (ft^2)
htail.lam = baseUAV.htail.lam          + (rand-0.5)*htail.sd_lam; % taper ratio of horizontal tail (btw 0 and 1 inclusive)
htail.lam_q = baseUAV.htail.lam_q      + (rand-0.5)*htail.sd_lam_q;  % horizontal tail sweep angle
% htail.lam_max = baseUAV.htail.lam_max  + (rand-0.5)*htail.sd_lam_max; % sweep of maximum thickness line 
% ^^ This value is derived parameter
htail.h = baseUAV.htail.h              + (rand-0.5)*htail.sd_h; % dist from head to 1/4 chord of horizontal tail (ft)
htail.Q   = baseUAV.htail.Q            + (rand-0.5)*htail.sd_Q; % horizontal tail interference factor 

% airfoil properties
airf_index_h      = randi(12); % randomly select airfoil
airfoilh.CL_alpha = airfoils(airf_index_h).CL_alpha_deg; % CL_alpha (/deg)
airfoilh.CLmax    = airfoils(airf_index_h).CLmax; % CLmax (wing)
airfoilh.maxthick = airfoils(airf_index_h).GEO.max_thick; % max thickness (relative to chord length)
airfoilh.maxthick_loc = airfoils(airf_index_h).GEO.max_thick_location; % max thickness location (relative to chord length)
airfoilh.perim    = airfoils(airf_index_h).GEO.perimeter; % circumference/perimeter of airfoil

htail.mtr = 0.24; % maximum thickness ratio
htail.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
htail.S_wet = 2.003*htail.S; % wet area for horizontal tail (ft^2)

% Vertical Tail -----------------------------------------------------------

% primay
vtail.S =  baseUAV.vtail.S         + (rand-0.5)*vtail.sd_S; % area (ft^2)
vtail.A = baseUAV.vtail.A          + (rand-0.5)*vtail.sd_A;   % aspect ratio (defined as b^2/S) 
vtail.lam = baseUAV.vtail.lam      + (rand-0.5)*vtail.sd_lam;  % taper ratio 
vtail.lam_q = baseUAV.vtail.lam_q  + (rand-0.5)*vtail.sd_lam_q;  % quarter chord sweep angle
vtail.lam_max = baseUAV.vtail.lam_max  + (rand-0.5)*vtail.sd_lam_max; % sweep of maximum thicknes line 
vtail.Q = baseUAV.vtail.Q          + (rand-0.5)*vtail.sd_Q;   % interference factor 
vtail.h = baseUAV.vtail.h          + (rand-0.5)*vtail.sd_h; % dist from head to 1/4 chord of vertical tail (ft)

% airfoil properties
airf_index_v      = randi(12); % randomly select airfoil
airfoilv.CL_alpha = airfoils(airf_index_v).CL_alpha_deg; % CL_alpha (/deg)
airfoilv.CLmax    = airfoils(airf_index_v).CLmax; % CLmax (wing)
airfoilv.maxthick = airfoils(airf_index_v).GEO.max_thick; % max thickness (relative to chord length)
airfoilv.maxthick_loc = airfoils(airf_index_v).GEO.max_thick_location; % max thickness location (relative to chord length)
airfoilv.perim    = airfoils(airf_index_v).GEO.perimeter; % circumference/perimeter of airfoil

vtail.mtr = 0.12; % maximum thickness ratio
vtail.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
vtail.S_wet = 2.003*vtail.S; % wet area for vertical tail (ft^2)

% Airfoil -----------------------------------------------------------------

%Randomly select airfoils from airfoil directory
airf_index       = randi(12); % randomly select airfoil
airfoil.CL_alpha = airfoils(airf_index).CL_alpha_deg; % CL_alpha (/deg)
airfoil.CLmax    = airfoils(airf_index).CLmax; % CLmax (wing)
airfoil.maxthick = airfoils(airf_index).GEO.max_thick; % max thickness (relative to chord length)
airfoil.maxthick_loc = airfoils(airf_index).GEO.max_thick_location; % max thickness location (relative to chord length)
airfoil.perim    = airfoils(airf_index).GEO.perimeter; % circumference/perimeter of airfoil

% Engine ------------------------------------------------------------------

%Randomly select engines from engine directory
eng_index = randi(9); % radomly select engine
engn.P_avail = engines(eng_index).P_avail; % load P_avail
engn.vol = engines(eng_index).vol; % load volume
engn.weight = engines(eng_index).weight; % load weight
engn.rpm = engines(eng_index).rpm; % load rpm

% Fuel System -------------------------------------------------------------

%   no items here

% Propeller ---------------------------------------------------------------

%   no items here

% Electronics/Payloads ----------------------------------------------------

payld.w_total = payld.w_total; % in lbs (no variation in payload)

% Surface Control ---------------------------------------------------------

%   no items here


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
