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
% Fuselage ----------------------------------------------------------------

fuse.L    = check_allowable(baseUAV.fuse.L,(rand-0.5)*fuse.sd_L,payld.length_TOTAL,100);    % fuselage max length
fuse.W    = check_allowable(baseUAV.fuse.W,(rand-0.5)*fuse.sd_W,13/12,5);    % fuselage max width 
fuse.D    = check_allowable(baseUAV.fuse.D,(rand-0.5)*fuse.sd_D,13/12,5);    % fuselage max depth

% Wing --------------------------------------------------------------------

wing.S     = check_allowable(baseUAV.wing.S,(rand-0.5)*wing.sd_S,0,30); % wing area (ft^2) 
wing.A     = check_allowable(baseUAV.wing.A,(rand-0.5)*wing.sd_A,0,30); % aspect ratio      
wing.lam   = check_allowable(baseUAV.wing.lam,(rand-0.5)*wing.sd_lam,0,1); % taper ratio (must be between 0 < 1)
wing.lam_q = check_allowable(baseUAV.wing.lam_q,(rand-0.5)*wing.sd_lam_q,0,1); % wing quarter chord sweep
wing.h     = check_allowable(baseUAV.wing.h,(rand-0.5)*wing.sd_h,0,fuse.L); % dist from head to wing 1/4 chord (ft) TODO: where in quarter chord?
wing.e = 0.8;     % Oswald's efficiency factor 
%wing.gamma = 2;   % wing dihedral (deg)  GLOBAL HAWK

% airfoil properties
airfoils = airfoil_to_wing(airfoils,wing.A, wing.e); % convert to wing props
airf_index_w       = randi(12); % randomly select airfoil
airfoilw.CL       = airfoils(airf_index_w).CL;  % Entire CL profile
airfoilw.Cd       = airfoils(airf_index_w).Cd;   % Drag due to airfoil choice
airfoilw.CLmax    = airfoils(airf_index_w).CLmax; % CLmax (wing)
airfoilw.maxthick = airfoils(airf_index_w).GEO.max_thick; % max thickness (relative to chord length)
airfoilw.maxthick_loc = airfoils(airf_index_w).GEO.max_thick_location; % max thickness location (relative to chord length)
airfoilw.perim    = airfoils(airf_index_w).GEO.perimeter; % circumference/perimeter of airfoil
airfoilw.a_w      = airfoils(airf_index_w).CL_alpha_rad; % lift curve slope of wing
airfoilw.alpha0   = airfoils(airf_index_w).alpha0; % zero angle of attack
airfoilw.alpha    = airfoils(airf_index_w).alpha; % angle of attack
airfoilw.Cm_ac    = airfoils(airf_index_w).CM;

wing.mtr = airfoilw.maxthick; % maximum thickness ratio
wing.mtl = airfoilw.maxthick_loc;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
wing.S_wet = 2.003*wing.S;  % wet area for wing (ft^2)

% Horizontal Tail ---------------------------------------------------------

% primary
htail.A = check_allowable(baseUAV.htail.A,(rand-0.5)*htail.sd_A,0,30);    % aspect ratio
htail.S = check_allowable(baseUAV.htail.S,(rand-0.5)*htail.sd_S,0,30); % area (ft^2)
htail.lam = check_allowable(baseUAV.htail.lam,(rand-0.5)*htail.sd_lam,0,1); % taper ratio of horizontal tail (btw 0 and 1 inclusive)
htail.lam_q = check_allowable(baseUAV.htail.lam_q,(rand-0.5)*htail.sd_lam_q,0,1);  % horizontal tail sweep angle
% htail.lam_max = baseUAV.htail.lam_max  + (rand-0.5)*htail.sd_lam_max; % sweep of maximum thickness line 
% ^^ This value is derived parameter
htail.h = check_allowable(baseUAV.htail.h,(rand-0.5)*htail.sd_h,1.05*wing.h,fuse.L); % dist from head to 1/4 chord of horizontal tail (ft)
htail.e = 0.8;      % Oswald Efficiency Factor

% airfoil properties
airfoils = airfoil_to_wing(airfoils,htail.A,htail.e);
airf_index_h      = randi(12); % randomly select airfoil
airfoilh.CL       = airfoils(airf_index_h).CL;  % Entire CL profile
airfoilh.Cd       = airfoils(airf_index_h).Cd;   % Drag due to airfoil choice
airfoilh.CLmax    = airfoils(airf_index_h).CLmax; % CLmax (wing)
airfoilh.maxthick = airfoils(airf_index_h).GEO.max_thick; % max thickness (relative to chord length)
airfoilh.maxthick_loc = airfoils(airf_index_h).GEO.max_thick_location; % max thickness location (relative to chord length)
airfoilh.perim    = airfoils(airf_index_h).GEO.perimeter; % circumference/perimeter of airfoil
airfoilh.a_t      = airfoils(airf_index_h).CL_alpha_rad; % lift curve slope of tail

htail.mtr = airfoilh.maxthick; % maximum thickness ratio
htail.mtl = airfoilh.maxthick_loc;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
htail.S_wet = 2.003*htail.S; % wet area for horizontal tail (ft^2)

% Vertical Tail -----------------------------------------------------------

% primay
vtail.S =  check_allowable(baseUAV.vtail.S,(rand-0.5)*vtail.sd_S,0,30); % area (ft^2)
vtail.A = check_allowable(baseUAV.vtail.A,(rand-0.5)*vtail.sd_A,0,30);   % aspect ratio (defined as b^2/S) 
vtail.lam = check_allowable(baseUAV.vtail.lam,(rand-0.5)*vtail.sd_lam,0,1);  % taper ratio 
vtail.lam_q = check_allowable(baseUAV.vtail.lam_q,(rand-0.5)*vtail.sd_lam_q,0,1);  % quarter chord sweep angle
vtail.h = check_allowable(baseUAV.vtail.h,(rand-0.5)*vtail.sd_h,1.05*wing.h,fuse.L); % dist from head to 1/4 chord of vertical tail (ft)
vtail.e = 0.8;  % Oswald Efficiency Factor

% airfoil properties
airfoils = airfoil_to_wing(airfoils,1.6*vtail.A,vtail.e);
airf_index_v      = randi(12); % randomly select airfoil
airfoilv.CL       = airfoils(airf_index_v).CL;  % Entire CL profile
airfoilv.Cd       = airfoils(airf_index_v).Cd;   % Drag due to airfoil choice
airfoilv.CL_alpha = airfoils(airf_index_v).CL_alpha_rad; % CL_alpha (/rad)
airfoilv.CLmax    = airfoils(airf_index_v).CLmax; % CLmax (wing)
airfoilv.maxthick = airfoils(airf_index_v).GEO.max_thick; % max thickness (relative to chord length)
airfoilv.maxthick_loc = airfoils(airf_index_v).GEO.max_thick_location; % max thickness location (relative to chord length)
airfoilv.perim    = airfoils(airf_index_v).GEO.perimeter; % circumference/perimeter of airfoil
airfoilv.a_v      = airfoils(airf_index_v).CL_alpha_rad; % lift curve slope of tail

vtail.mtr = airfoilv.maxthick; % maximum thickness ratio
vtail.mtl = airfoilv.maxthick_loc;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
vtail.S_wet = 2.003*vtail.S; % wet area for vertical tail (ft^2)


% Engine ------------------------------------------------------------------

%Randomly select engines from engine directory
eng_index = randi(9); % radomly select engine
engn.P_avail = engines(eng_index).P_avail; % load P_avail
engn.vol = engines(eng_index).vol; % load volume
engn.weight = engines(eng_index).weight; % load weight
engn.rpm = engines(eng_index).rpm; % load rpm

% Fuel --------------------------------------------------------------------

fuel.cp = 0.85/550/3600;  %[1/ft] specific fuel consumption
fuel.rho = 6.073; %[lbm/gallon] density of fuel for octane gas
fuel.W = check_allowable(baseUAV.fuel.W,(rand-0.5)*fuel.sd_W,0,100);
fuel.V = fuel.W/fuel.rho; %[gallon] volume of fuel
fuel.V = fuel.V*gallon2ft3; %[ft^3] volume of fuel

% Fuel System -------------------------------------------------------------

fsys.W = (fuel.V-0.00424)*1726.635+65; % in grams (emperically derived)
fsys.W = fsys.W*g2lb; % in pounds
fsys.length = fsys.W/(pi*(fuse.D/2)^2); % fuel tank length

% Propeller ---------------------------------------------------------------

%   no items here

% Electronics/Payloads ----------------------------------------------------

payld.w_total = payld.w_total; % in lbs (no variation in payload)

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
wing.x_cg = wing.h;
wing.z_cg = baseUAV.wing.z_cg;

wing.h_q_t = wing.h+tand(wing.lam_q)*wing.b/2;  % wing c/4 tip location (ft)
wing.h_l_r = wing.h-wing.c_r/4;                 % wing leading edge root location (ft)
wing.h_l_t = wing.h_q_t-wing.c_t/4;             % wing leading edge tip location (ft)
wing.lam_l = atand((wing.h_l_t-wing.h_l_r)/(wing.b/2)); % wing leading edge sweep location (deg)

% Fuselage ----------------------------------------------------------------

fuse.x_cg = fuse.L*0.5; % Assume 50% fuselage length
fuse.z_cg = baseUAV.fuse.z_cg;
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
htail.t_t = htail.c_t*htail.mtr;   % max thickness ratio (ft)s
htail.l_T = htail.h-wing.h;      % distance from wing 1/4 MAC to tail 1/4 MAC (ft)
htail.x_cg = htail.h;
htail.z_cg = baseUAV.htail.z_cg;

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
vtail.x_cg = vtail.h;
vtail.z_cg = baseUAV.vtail.z_cg;

vtail.h_q_t = vtail.h+tand(vtail.lam_q)*vtail.b;
vtail.h_l_r = vtail.h    -vtail.c_r/4;
vtail.h_l_t = vtail.h_q_t-vtail.c_t/4;
vtail.lam_l = atand((vtail.h_l_t-vtail.h_l_r)/(vtail.b/2));

% Engine ------------------------------------------------------------------

engn.x_cg = 0.95*fuse.L; % engine CG location (assume 95% of fuselage)
engn.z_cg = baseUAV.engn.z_cg;

% Fuel System -------------------------------------------------------------

fsys.x_cg = check_allowable(baseUAV.fsys.x_cg,(rand-0.5)*fsys.sd_x_cg,0,0.9*fuse.L) ; % fuel system CG location (ft)
fsys.z_cg = baseUAV.fsys.z_cg;

% Propeller ---------------------------------------------------------------

prop.h    = fuse.L; % beginning of prop is located at end of fuselge
prop.x_cg = (1+0.025)*fuse.L; % propeller CG location (ft)
prop.z_cg = baseUAV.prop.z_cg;

% Electronics/Payloads ----------------------------------------------------

% CG locations
payld.x_cg_EOIR  = 0.025*fuse.L; % EO/IR CG location (ft)
payld.x_cg_SAR   = 0.05*fuse.L;  % Synthetic Apperature Radar CG location (ft)
payld.x_cg_LiDAR = check_allowable(0.075*fuse.L,(rand-0.5)*payld.sd_x_cg_LiDAR,0,0.9*fuse.L); % LiDAR CG location (ft)
payld.x_cg_ANT   = check_allowable(0.05*fuse.L,(rand-0.5)*payld.sd_x_cg_ANT,0,0.9*fuse.L);  % UHF/VHF antenna location (ft)
payld.x_cg_WR    = check_allowable(0.05*fuse.L,(rand-0.5)*payld.sd_x_cg_WR,0,0.9*fuse.L);  % WaveRelay CG location (ft)
payld.x_cg_IMU   = check_allowable(0.05*fuse.L,(rand-0.5)*payld.sd_x_cg_IMU,0,0.9*fuse.L);  % IMU CG location (ft)

payld.z_cg_EOIR  = baseUAV.payld.z_cg_EOIR;
payld.z_cg_SAR   = baseUAV.payld.z_cg_SAR;
payld.z_cg_LiDAR = baseUAV.payld.z_cg_LiDAR;
payld.z_cg_ANT   = baseUAV.payld.z_cg_ANT;
payld.z_cg_WR    = baseUAV.payld.z_cg_WR;
payld.z_cg_IMU   = baseUAV.payld.z_cg_IMU;

% ============================ CONTROL SURFACE ============================
% Aileron -----------------------------------------------------------------

% primary
sfcl.ail.S    = check_allowable(baseUAV.sfcl.ail.S,(rand-0.5)*sfcl.ail.sd_S,0.03*wing.S,0.12*wing.S); % aileron area
sfcl.ail.c    = check_allowable(baseUAV.sfcl.ail.c,(rand-0.5)*sfcl.ail.sd_c,0.15*wing.c,0.3*wing.c); % aileron chord
sfcl.ail.sc_W = 0.09193276; % weight of servo in lbs

% derived
sfcl.ail.b    = sfcl.ail.S/sfcl.ail.c; % aileron length

% Rudder ------------------------------------------------------------------

% primary
sfcl.rudd.S    = check_allowable(baseUAV.sfcl.rudd.S,(rand-0.5)*sfcl.rudd.sd_S,0.15*vtail.S,0.3*vtail.S); % rudder area
sfcl.rudd.c    = check_allowable(baseUAV.sfcl.rudd.c,(rand-0.5)*sfcl.rudd.sd_c,0.15*vtail.c,0.4*vtail.c); % rudder chord
sfcl.rudd.sc_W = 0.0198416; % weight of servo in lbs

% derived
sfcl.rudd.b    = sfcl.rudd.S/sfcl.rudd.c; % rudder length

% Elevator ----------------------------------------------------------------

% primary
sfcl.elev.S    = check_allowable(baseUAV.sfcl.elev.S,(rand-0.5)*sfcl.elev.sd_S,0.15*htail.S,0.4*htail.S); % elevator area
sfcl.elev.c    = check_allowable(baseUAV.sfcl.elev.c,(rand-0.5)*sfcl.elev.sd_c,0.2*htail.c,0.4*htail.c); % elevator chord
sfcl.elev.sc_W = 0.0198416; % weight of servo in lbs

% derived
sfcl.elev.b    = sfcl.elev.S/sfcl.elev.c; % elevator length

sfcl.ail.x_cg  = wing.x_cg+wing.c*(0.85-0.25); % assume 85% of average wing chord
sfcl.elev.x_cg = htail.x_cg+htail.c*(0.85-0.25); % assume 85% of average wing chord
sfcl.rudd.x_cg = vtail.x_cg+htail.c*(0.85-0.25); % assume 85% of average wing chord

sfcl.ail.z_cg  = wing.z_cg;
sfcl.elev.z_cg = htail.z_cg;
sfcl.rudd.z_cg = vtail.z_cg;
