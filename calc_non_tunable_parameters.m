% calc_non_tunable_parameters

% Wing --------------------------------------------------------------------

wing.e = 0.8;     % Oswald's efficiency factor 

% airfoil properties for wing
airfoils = airfoil_to_wing(airfoils,wing.A, wing.e); % convert to wing props
airfoilw.CL       = airfoils(airfoilw.ind).CL;  % Entire CL profile
airfoilw.Cd       = airfoils(airfoilw.ind).Cd;   % Drag due to airfoil choice
airfoilw.CLmax    = airfoils(airfoilw.ind).CLmax; % CLmax (wing)
airfoilw.maxthick = airfoils(airfoilw.ind).GEO.max_thick; % max thickness (relative to chord length)
airfoilw.maxthick_loc = airfoils(airfoilw.ind).GEO.max_thick_location; % max thickness location (relative to chord length)
airfoilw.perim    = airfoils(airfoilw.ind).GEO.perimeter; % circumference/perimeter of airfoil
airfoilw.a_w      = airfoils(airfoilw.ind).CL_alpha_rad; % lift curve slope of wing
airfoilw.alpha0   = airfoils(airfoilw.ind).alpha0; % zero angle of attack
airfoilw.alpha    = airfoils(airfoilw.ind).alpha; % angle of attack
airfoilw.Cm_ac    = airfoils(airfoilw.ind).CM;

wing.mtr = airfoilw.maxthick; % maximum thickness ratio
wing.mtl = airfoilw.maxthick_loc;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
wing.S_wet = 2.1*wing.S;  % wet area for wing (ft^2)

wing.b = sqrt(wing.A*wing.S);     % wing span (ft)
wing.c = wing.S/wing.b;           % average wing cord (ft)   TODO: change name to wing.c_ave
wing.c_r = 2*wing.c/(1+wing.lam); % wing root chord length (ft)
wing.c_t = wing.c_r*wing.lam;     % wing tip chord length (ft)
wing.t_r = wing.c_r*wing.mtr;     % max thickness at root (ft)
wing.t_t = wing.c_t*wing.mtr;     % max thickness at tip (ft)
wing.t = wing.c*wing.mtr;         % average thickness of wing (ft) TODO: change name to wing.t_ave
wing.K_i = 1/(pi*wing.A*wing.e);  % Drag coefficient
wing.x_cg = wing.h;
wing.z_cg = (0.5*fuse.D)+(0.5*wing.t);

wing.h_q_t = wing.h+tand(wing.lam_q)*wing.b/2;  % wing c/4 tip location (ft)
wing.h_l_r = wing.h-wing.c_r/4;                 % wing leading edge root location (ft)
wing.h_l_t = wing.h_q_t-wing.c_t/4;             % wing leading edge tip location (ft)
wing.lam_l = atand((wing.h_l_t-wing.h_l_r)/(wing.b/2)); % wing leading edge sweep location (deg)

% Fuselage ----------------------------------------------------------------

fuse.r = fuse.L/fuse.W; %fuselage ratio

% arbitrarily defined 
fuse.S_wet = pi*0.5*((fuse.L*fuse.W)+(fuse.L*fuse.W));%wet fuselage area (ft^2)

% Horizontal Tail ---------------------------------------------------------

htail.e = 0.8;      % Oswald Efficiency Factor

% airfoil properties for horizontal tail
airfoils = airfoil_to_wing(airfoils,htail.A,htail.e);

airfoilh.CL       = airfoils(aifoilh.ind).CL;  % Entire CL profile
airfoilh.Cd       = airfoils(aifoilh.ind).Cd;   % Drag due to airfoil choice
airfoilh.CLmax    = airfoils(aifoilh.ind).CLmax; % CLmax (wing)
airfoilh.maxthick = airfoils(aifoilh.ind).GEO.max_thick; % max thickness (relative to chord length)
airfoilh.maxthick_loc = airfoils(aifoilh.ind).GEO.max_thick_location; % max thickness location (relative to chord length)
airfoilh.perim    = airfoils(aifoilh.ind).GEO.perimeter; % circumference/perimeter of airfoil
airfoilh.a_t      = airfoils(aifoilh.ind).CL_alpha_rad; % lift curve slope of tail

htail.mtr = airfoilh.maxthick; % maximum thickness ratio
htail.mtl = airfoilh.maxthick_loc;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
htail.S_wet = 2.003*htail.S; % wet area for horizontal tail (ft^2)

htail.b = sqrt(htail.A*htail.S); % span (ft)
htail.c = htail.S/htail.b;       % average chord length (ft) TODO: change name to htail.c_ave
htail.c_r = 2*htail.c/(1+htail.lam); % horizontal tail root chord length (ft)
htail.c_t = htail.c_r*htail.lam;     % wing tip chord length (ft)
% htail.t   = htail.c*htail.mtr;     % max thickness (average) (ft) TODO: change name to htail.t_ave
htail.t_r = htail.c_r*htail.mtr;   % max thickness at root (ft)
htail.t_t = htail.c_t*htail.mtr;   % max thickness ratio (ft)s
htail.l_T = htail.h-wing.h;      % distance from wing 1/4 MAC to tail 1/4 MAC (ft)

htail.h_q_t = htail.h+tand(htail.lam_q)*htail.b/2;
htail.h_l_r = htail.h    -htail.c_r/4;
htail.h_l_t = htail.h_q_t-htail.c_t/4;
htail.lam_l = atand((htail.h_l_t-htail.h_l_r)/(htail.b/2));


% Vertical Tail -----------------------------------------------------------

vtail.e = 0.8;  % Oswald Efficiency Factor

% airfoil properties
airfoils = airfoil_to_wing(airfoils,1.6*vtail.A,vtail.e);
airfoilv.CL       = airfoils(airfoilv.ind).CL;  % Entire CL profile
airfoilv.Cd       = airfoils(airfoilv.ind).Cd;   % Drag due to airfoil choice
airfoilv.CL_alpha = airfoils(airfoilv.ind).CL_alpha_rad; % CL_alpha (/rad)
airfoilv.CLmax    = airfoils(airfoilv.ind).CLmax; % CLmax (wing)
airfoilv.maxthick = airfoils(airfoilv.ind).GEO.max_thick; % max thickness (relative to chord length)
airfoilv.maxthick_loc = airfoils(airfoilv.ind).GEO.max_thick_location; % max thickness location (relative to chord length)
airfoilv.perim    = airfoils(airfoilv.ind).GEO.perimeter; % circumference/perimeter of airfoil
airfoilv.a_v      = airfoils(airfoilv.ind).CL_alpha_rad; % lift curve slope of tail

vtail.mtr = airfoilv.maxthick; % maximum thickness ratio
vtail.mtl = airfoilv.maxthick_loc;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
vtail.S_wet = 2.003*vtail.S; % wet area for vertical tail (ft^2)

vtail.b   = sqrt(vtail.A*vtail.S); % span (ft)
vtail.c   = vtail.S/vtail.b;       % average chord length (ft)
vtail.c_r = 2*vtail.c/(1+vtail.lam); % horizontal tail root chord length (ft)
vtail.c_t = vtail.c_r*vtail.lam;     % wing tip chord length (ft)
vtail.t = vtail.c*vtail.mtr;     % max root thickness average(ft)
vtail.t_r =vtail.c_r*vtail.mtr;  % max root thickness at root (ft)
vtail.t_t = vtail.c_t*vtail.mtr; % max root thickness at tip (ft)

vtail.h_q_t = vtail.h+tand(vtail.lam_q)*vtail.b;
vtail.h_l_r = vtail.h    -vtail.c_r/4;
vtail.h_l_t = vtail.h_q_t-vtail.c_t/4;
vtail.lam_l = atand((vtail.h_l_t-vtail.h_l_r)/(vtail.b/2));

% Engine ------------------------------------------------------------------

engn.N_E    = 1;         % number of engines

engn.length  = engines(engn.ind).l;  % cylinder part of the engine length (ft)
engn.HP      = engines(engn.ind).P_avail; % power available (hp) (p_max*engine_efficiency)
engn.vol     = engines(engn.ind).vol;     % engine total volume (ft^3)
engn.W_bare  = engines(engn.ind).weight; % load weight
engn.rpm     = engines(engn.ind).rpm; % load rpm

% Fuel --------------------------------------------------------------------

fuel.cp = 0.85/550/3600;  %[1/ft] specific fuel consumption
fuel.rho = 6.073; %[lbm/gallon] density of fuel for octane gas

fuel.V = fuel.W/fuel.rho; %[gallon] volume of fuel
fuel.V = fuel.V*gallon2ft3; %[ft^3] volume of fuel

% Propeller ---------------------------------------------------------------

prop.h    = fuse.L; % beginning of prop is located at end of fuselge

% Fuel System -------------------------------------------------------------

fsys.int = 0.00; %percent of fuel tanks that are integral TODO: find out how much is integral
fsys.N_T = 1;    %number of separate fuel tanks

fsys.W = (fuel.V-0.00424)*1726.635+65; % in grams (emperically derived)
fsys.W = fsys.W*g2lb; % in pounds
fsys.length = fuel.V/(pi*(fuse.D-(1/12)^2));  % ft

% Electronics/Payloads ----------------------------------------------------

payld.w_EOIR  = 12.57; % EO/IR weight(lbs)
payld.w_SAR   = 2;     % Synthetic Apperature Radar weight(lbs)
payld.w_LiDAR = 1;     % LiDAR weight(lbs)
payld.w_ANT   = 0.15;  % UHF/VHF antenna weight (lbs)
payld.w_WR    = 0.25;  % WaveRelay weight (lbs)
payld.w_IMU   = 0.11;  % Ellipse E mini (lbs)

% Pyaload/Avionics total weight
payld.w_total = payld.w_EOIR + payld.w_SAR + payld.w_LiDAR + payld.w_ANT ...
                + payld.w_WR + payld.w_IMU;

% Aileron ------------------------------------------------------------------

sfcl.ail.sc_W = 0.09193276; % weight of servo in lbs
sfcl.ail.b    = sfcl.ail.S/sfcl.ail.c; % aileron length

% Rudder ------------------------------------------------------------------

sfcl.rudd.sc_W = 0.0198416; % weight of servo in lbs
sfcl.rudd.b    = sfcl.rudd.S/sfcl.rudd.c; % rudder length

% Elevator ----------------------------------------------------------------

sfcl.elev.sc_W = 0.0198416; % weight of servo in lbs
sfcl.elev.b    = sfcl.elev.S/sfcl.elev.c; % elevator length

% CG locations ------------------------------------------------------------

wing.x_cg  = wing.h;   % (ft)
fuse.x_cg  = fuse.L/2; % (ft)
htail.x_cg = htail.h;  % (ft)
vtail.x_cg = vtail.h;  % (ft)
engn.x_cg = 0.95*fuse.L;  % (ft)
fsys.x_cg = 0.65*fuse.L;  % (ft)
% fsys.x_cg = check_allowable(baseUAV.fsys.x_cg,(rand-0.5)*fsys.sd_x_cg,0,0.9*fuse.L) ; % fuel system CG location (ft)
prop.x_cg = 1.025*fuse.L; % (ft) TODO: why 2.5%?
payld.x_cg_EOIR  = 0.025*fuse.L; % (ft)
payld.x_cg_SAR   = 0.05*fuse.L;  % (ft)
sfcl.ail.x_cg  = wing.x_cg+wing.c*(0.85-0.25); % assume 85% of average wing chord
sfcl.elev.x_cg = htail.x_cg+htail.c*(0.85-0.25); % assume 85% of average wing chord
sfcl.rudd.x_cg = vtail.x_cg+htail.c*(0.85-0.25); % assume 85% of average wing chord


wing.z_cg = (0.5*fuse.D)+(0.5*wing.t); % Assume z_cg is middle of average thickness (ft)
fuse.z_cg = fuse.D/2;
htail.z_cg = (0.5*fuse.D)+(0.5*wing.t); % tail z-CG location (ft)
vtail.z_cg = fuse.D+(0.5*vtail.c);
engn.z_cg = 0.5*fuse.D;
fsys.z_cg = 0.5*fuse.D;
prop.z_cg = engn.z_cg;
payld.z_cg_EOIR = 0.5*fuse.D;
payld.z_cg_SAR  = 0.5*fuse.D;
payld.z_cg_LiDAR = 0.5*fuse.D;
payld.z_cg_ANT = 0;
payld.z_cg_WR = 0.5*fuse.D;
payld.z_cg_IMU = 0.5*fuse.D;
sfcl.ail.z_cg  = wing.z_cg;
sfcl.elev.z_cg = htail.z_cg;
sfcl.rudd.z_cg = vtail.z_cg;

% Payload lengths x-direction
payld.length_EOIR   = 10.04*in2ft;  % (ft)
payld.length_SAR    = 6.2*in2ft;    % (ft)
payld.length_LiDAR  = 6.47*in2ft;   % (ft)
payld.length_WR_IMU = 3.9*in2ft;    % (ft)

% Payload-Unit
payld.length_TOTAL = (0.2*m2ft) + payld.length_EOIR + payld.length_SAR + ...
                    payld.length_LiDAR; % (ft)