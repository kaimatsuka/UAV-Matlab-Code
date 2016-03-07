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

VARY_WING  = 1;
VARY_FUSE  = 1;
VARY_HTAIL = 1;
VARY_VTAIL = 1;
VARY_FUEL  = 1;
VARY_ENGN  = 1;
VARY_PROP  = 0;
VARY_SFCL  = 0;
VARY_CG    = 0;

DO_NOT_VARY_ANYTHING = 0;

if DO_NOT_VARY_ANYTHING
    VARY_WING  = 0;
    VARY_FUSE  = 0;
    VARY_HTAIL = 0;
    VARY_VTAIL = 0;
    VARY_FUEL  = 0;
    VARY_ENGN  = 0;
    VARY_PROP  = 0;
    VARY_SFCL  = 0;
    VARY_CG    = 0;
end


% Wing --------------------------------------------------------------------

if VARY_WING
    wing.S     = check_allowable(baseUAV.wing.S,(rand-0.5)*wing.sd_S,0,30); % wing area (ft^2) 
    wing.A     = check_allowable(baseUAV.wing.A,(rand-0.5)*wing.sd_A,0,30); % aspect ratio      
    wing.lam   = check_allowable(baseUAV.wing.lam,(rand-0.5)*wing.sd_lam,0,1); % taper ratio (must be between 0 < 1)
    wing.lam_q = check_allowable(baseUAV.wing.lam_q,(rand-0.5)*wing.sd_lam_q,0,1); % wing quarter chord sweep
    wing.h_q   = check_allowable(baseUAV.wing.h_q,(rand-0.5)*wing.sd_h,0,fuse.L); % dist from head to wing 1/4 chord (ft) TODO: where in quarter chord?

    % randomly pick airfoil
    airfoilw.ind  = randi(12); % randomly select airfoil
else
    wing.S     = baseUAV.wing.S;
    wing.A     = baseUAV.wing.A;
    wing.lam   = baseUAV.wing.lam;
    wing.lam_q = baseUAV.wing.lam_q;
    wing.h_q   = baseUAV.wing.h_q;     % dist from head to wing 1/4 chord at root (ft)

    % randomly pick airfoil
    airfoilw.ind  = baseUAV.airfoilw.ind; % by default, select airfoil 2
end

% Fuselage ----------------------------------------------------------------

if VARY_FUSE
    fuse.L    = check_allowable(baseUAV.fuse.L,(rand-0.5)*fuse.sd_L,payld.length_TOTAL,100);    % fuselage max length
    fuse.W    = check_allowable(baseUAV.fuse.W,(rand-0.5)*fuse.sd_W,13/12,5);    % fuselage max width 
    fuse.D    = check_allowable(baseUAV.fuse.D,(rand-0.5)*fuse.sd_D,13/12,5);    % fuselage max depth
else
    fuse.L    = baseUAV.fuse.L;
    fuse.W    = baseUAV.fuse.W;
    fuse.D    = baseUAV.fuse.D;
end
    
% Horizontal Tail ---------------------------------------------------------

if VARY_HTAIL
    % primary
    htail.A = check_allowable(baseUAV.htail.A,(rand-0.5)*htail.sd_A,0,30);    % aspect ratio
    htail.S = check_allowable(baseUAV.htail.S,(rand-0.5)*htail.sd_S,0,30); % area (ft^2)
    htail.lam = check_allowable(baseUAV.htail.lam,(rand-0.5)*htail.sd_lam,0,1); % taper ratio of horizontal tail (btw 0 and 1 inclusive)
    htail.lam_q = check_allowable(baseUAV.htail.lam_q,(rand-0.5)*htail.sd_lam_q,0,1);  % horizontal tail sweep angle
    htail.h = check_allowable(baseUAV.htail.h,(rand-0.5)*htail.sd_h,1.05*wing.h_q,fuse.L); % dist from head to 1/4 chord of horizontal tail (ft)
    % htail.lam_max = baseUAV.htail.lam_max  + (rand-0.5)*htail.sd_lam_max; % sweep of maximum thickness line 
    % ^^ This value is derived parameter

    % randomly pick airfoil
    airfoilh.ind = randi(12); % randomly select airfoil
else
    % primary
    htail.A     = baseUAV.htail.A;
    htail.S     = baseUAV.htail.S;
    htail.lam   = baseUAV.htail.lam;
    htail.lam_q = baseUAV.htail.lam_q;
    htail.h     = baseUAV.htail.h;
    % htail.lam_max = baseUAV.htail.lam_max  + (rand-0.5)*htail.sd_lam_max; % sweep of maximum thickness line 
    % ^^ This value is derived parameter

    % randomly pick airfoil
    airfoilh.ind = baseUAV.airfoilh.ind; % randomly select airfoil
end

% Vertical Tail -----------------------------------------------------------

if VARY_VTAIL
    % primay
    vtail.S =  check_allowable(baseUAV.vtail.S,(rand-0.5)*vtail.sd_S,0,30); % area (ft^2)
    vtail.A = check_allowable(baseUAV.vtail.A,(rand-0.5)*vtail.sd_A,0,30);   % aspect ratio (defined as b^2/S) 
    vtail.lam = check_allowable(baseUAV.vtail.lam,(rand-0.5)*vtail.sd_lam,0,1);  % taper ratio 
    vtail.lam_q = check_allowable(baseUAV.vtail.lam_q,(rand-0.5)*vtail.sd_lam_q,0,1);  % quarter chord sweep angle
    vtail.h = check_allowable(baseUAV.vtail.h,(rand-0.5)*vtail.sd_h,1.05*wing.h_q,fuse.L); % dist from head to 1/4 chord of vertical tail (ft)

    % randomly pick airfoil
    airfoilv.ind = randi(12); % randomly select airfoil
else
    % primay
    vtail.S = baseUAV.vtail.S;
    vtail.A = baseUAV.vtail.A;
    vtail.lam = baseUAV.vtail.lam;
    vtail.lam_q = baseUAV.vtail.lam_q;
    vtail.h = baseUAV.vtail.h;

    % randomly pick airfoil
    airfoilv.ind = baseUAV.airfoilv.ind; % randomly select airfoil
end

% Fuel --------------------------------------------------------------------

if VARY_FUEL
    fuel.W = check_allowable(baseUAV.fuel.W,(rand-0.5)*fuel.sd_W,0,9999); %[lb] fuel weight (calculated using test_fuel.m file)
else
    fuel.W = baseUAV.fuel.W;
end

% Fuel System -------------------------------------------------------------

%   no items here

% Engine ------------------------------------------------------------------

if VARY_ENGN
    %Randomly select engines from engine directory
    engn.ind = randi(9); % index of engines
else
    engn.ind = 6;
end

% Propeller ---------------------------------------------------------------

if VARY_PROP
    prop.D     = check_allowable(baseUAV.prop.D,(rand-0.5)*prop.sd_D,0,9999); %[ft]
    prop.pitch = check_allowable(baseUAV.prop.pitch,(rand-0.5)*prop.sd_pitch,0,9999); %[ft]
else
    prop.D     = baseUAV.prop.D;
    prop.pitch = baseUAV.prop.pitch;
end

% Electronics/Payloads ----------------------------------------------------


if VARY_SFCL
    sfcl.ail.S    = check_allowable(baseUAV.sfcl.ail.S,(rand-0.5)*sfcl.ail.sd_S,0.03*wing.S,0.12*wing.S); % aileron area
    sfcl.ail.c    = check_allowable(baseUAV.sfcl.ail.c,(rand-0.5)*sfcl.ail.sd_c,0.15*wing.c,0.3*wing.c); % aileron chord

    sfcl.rudd.S    = check_allowable(baseUAV.sfcl.rudd.S,(rand-0.5)*sfcl.rudd.sd_S,0.15*vtail.S,0.3*vtail.S); % rudder area
    sfcl.rudd.c    = check_allowable(baseUAV.sfcl.rudd.c,(rand-0.5)*sfcl.rudd.sd_c,0.15*vtail.c,0.4*vtail.c); % rudder chord

    sfcl.elev.S    = check_allowable(baseUAV.sfcl.elev.S,(rand-0.5)*sfcl.elev.sd_S,0.15*htail.S,0.4*htail.S); % elevator area
    sfcl.elev.c    = check_allowable(baseUAV.sfcl.elev.c,(rand-0.5)*sfcl.elev.sd_c,0.2*htail.c,0.4*htail.c); % elevator chord
else
    
    sfcl.ail.S  = baseUAV.sfcl.ail.S;
    sfcl.ail.c  = baseUAV.sfcl.ail.c;

    sfcl.rudd.S = baseUAV.sfcl.rudd.S;
    sfcl.rudd.c = baseUAV.sfcl.rudd.c;

    sfcl.elev.S = baseUAV.sfcl.elev.S;
    sfcl.elev.c = baseUAV.sfcl.elev.c;
end


if VARY_CG
    payld.x_cg_LiDAR = check_allowable(baseUAV.payld.x_cg_LiDAR,(rand-0.5)*payld.sd_x_cg_LiDAR,0,0.9*fuse.L); % LiDAR CG location (ft)
    payld.x_cg_ANT   = check_allowable(baseUAV.payld.x_cg_ANT,(rand-0.5)*payld.sd_x_cg_ANT,0,0.9*fuse.L);  % UHF/VHF antenna location (ft)
    payld.x_cg_WR    = check_allowable(baseUAV.payld.x_cg_WR,(rand-0.5)*payld.sd_x_cg_WR,0,0.9*fuse.L);  % WaveRelay CG location (ft)
    payld.x_cg_IMU   = check_allowable(baseUAV.payld.x_cg_WR,(rand-0.5)*payld.sd_x_cg_IMU,0,0.9*fuse.L);  % IMU CG location (ft)
else
    payld.x_cg_LiDAR = baseUAV.payld.x_cg_LiDAR;
    payld.x_cg_ANT   = baseUAV.payld.x_cg_ANT;
    payld.x_cg_WR    = baseUAV.payld.x_cg_WR;
    payld.x_cg_IMU   = baseUAV.payld.x_cg_IMU;
end

%% Calculate derived geometries/variables

calc_non_tunable_parameters
