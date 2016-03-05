% load_variation_parameters.m
%
% DESCRIPTION:
%   This file assigns varation parameters (standard deviation) associated
%   with each of adjustable UAV parameters.
%
% INPUT:
%   N/A
%
% OUTPUT:
%   variation parameters
%
% REVISION HISTORY:
%   02/10: File created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set 1 to vary, set 0 to not
VARY_WING  = 1;
VARY_FUSE  = 1;
VARY_HTAIL = 1;
VARY_VTAIL = 1;
VARY_SFCL  = 1;
VARY_PROP  = 1;
VARY_FUEL  = 1;
VARY_CG    = 1;

% Wing --------------------------------------------------------------------
if VARY_WING
    wing.sd_A     = 10;  % wing aspect ratio variation std dev 
    wing.sd_S     = 15;   % wing area variation (ft^2)
    wing.sd_lam   = 0.5; % wing taper ratio variation(n/a)
    wing.sd_lam_q = 0;   % wing quarter chord sweep variation (deg)
    wing.sd_h     = 0.1;   % head to 1/4 chord distance variation (ft)
else
    wing.sd_A     = 0;  % wing aspect ratio variation std dev 
    wing.sd_S     = 0;   % wing area variation (ft^2)
    wing.sd_lam   = 0; % wing taper ratio variation(n/a)
    wing.sd_lam_q = 0;   % wing quarter chord sweep variation (deg)
    wing.sd_h     = 0;   % head to 1/4 chord distance variation (ft)
end

% Fuselage ----------------------------------------------------------------

if VARY_FUSE
    fuse.sd_L     = 8; % fuselage max length variation (ft)
    fuse.sd_W     = 0; % fuselage max width variation (ft)
    fuse.sd_D     = 0; % fuselage max depth variation (ft)
else
    fuse.sd_L     = 0; % fuselage max length variation (ft)
    fuse.sd_W     = 0; % fuselage max width variation (ft)
    fuse.sd_D     = 0; % fuselage max depth variation (ft)
end

% Horizontal Tail ---------------------------------------------------------

if VARY_HTAIL
    %primary
    htail.sd_S       = 0; % htail area variation (ft^2)
    htail.sd_A       = 5; % htail aspect ratio variation std dev
    htail.sd_lam     = 0; % htail taper ratio variation
    htail.sd_lam_q   = 0; % htail quarter chord sweep variation 
    % htail.sd_lam_max = 0; % htail sweep of maximum thickness line variation
    htail.sd_h       = 0.15; % htail dist from head to 1/4 chord of variation (ft)
else
    %primary
    htail.sd_S       = 0; % htail area variation (ft^2)
    htail.sd_A       = 0; % htail aspect ratio variation std dev
    htail.sd_lam     = 0; % htail taper ratio variation
    htail.sd_lam_q   = 0; % htail quarter chord sweep variation 
    % htail.sd_lam_max = 0; % htail sweep of maximum thickness line variation
    htail.sd_h       = 0; % htail dist from head to 1/4 chord of variation (ft)
end

% Vertical Tail -----------------------------------------------------------

if VARY_VTAIL
    vtail.sd_S       = 10; % vtail area variation (ft^2)
    vtail.sd_A       = 5; % vtail aspect ratio variation 
    vtail.sd_lam     = 0; % vtail taper ratio variation
    vtail.sd_lam_q   = 0; % vtail quarter chord sweep variation
    % vtail.sd_lam_max = 0; % vtail sweep of maximum thickness line variation
    vtail.sd_h       = 0.15; % vtail htail interference factor variation
else
    vtail.sd_S       = 0; % vtail area variation (ft^2)
    vtail.sd_A       = 0; % vtail aspect ratio variation 
    vtail.sd_lam     = 0; % vtail taper ratio variation
    vtail.sd_lam_q   = 0; % vtail quarter chord sweep variation
    % vtail.sd_lam_max = 0; % vtail sweep of maximum thickness line variation
    vtail.sd_h       = 0; % vtail htail interference factor variation
end

% Surface Controls --------------------------------------------------------

if VARY_SFCL 

    sfcl.ail.sd_S   = 0.5;
    sfcl.ail.sd_c   = 0.05;

    sfcl.rudd.sd_S  = 0.25;
    sfcl.rudd.sd_c  = 0.05;

    sfcl.elev.sd_S  = 0.25;
    sfcl.elev.sd_c  = 0.05;
else
    sfcl.ail.sd_S   = 0;
    sfcl.ail.sd_c   = 0;
    sfcl.rudd.sd_S  = 0;
    sfcl.rudd.sd_c  = 0;

    sfcl.elev.sd_S  = 0;
    sfcl.elev.sd_c  = 0;
end


% Propeller ---------------------------------------------------------------

if VARY_PROP 
    prop.sd_D     = 1.6; % propeller diameter variation (ft)
    prop.sd_pitch = 20;  % propeller pitch (ft)
else
    prop.sd_D     = 0;    % propeller diameter variation (ft)
    prop.sd_pitch = 0; % propeller pitch (ft)
end

% Fuel --------------------------------------------------------------------

if VARY_FUEL 
    fuel.sd_W = 2;   % lbs
else
    fuel.sd_W = 0;   % lbs
end

% CG locations ------------------------------------------------------------

if VARY_CG
    wing.sd_x_cg       = 4;  % (ft)
    fuse.sd_x_cg       = 4;  % (ft)
    htail.sd_x_cg      = 4; % (ft)
    vtail.sd_x_cg      = 4; % (ft)
    engn.sd_x_cg       = 0;  % (ft)
    fsys.sd_x_cg       = 4;  % (ft)
    prop.sd_x_cg       = 0;  % (ft)
    payld.sd_x_cg_EOIR = 0;  % (ft)
    payld.sd_x_cg_SAR  = 0;   % (ft)
    payld.sd_x_cg_LiDAR = 0; % (ft)
    payld.sd_x_cg_ANT  = 0;   % (ft)
    payld.sd_x_cg_WR   = 0;   % (ft)
    payld.sd_x_cg_IMU  = 0;   % (ft)
    sfcl.sd_x_cg_wing  = 0;  % (ft)
    sfcl.sd_x_cg_htail = 0;  % (ft)
    sfcl.sd_x_cg_vtail = 0;  % (ft)
else
    wing.sd_x_cg       = 0;  % (ft)
    fuse.sd_x_cg       = 0;  % (ft)
    htail.sd_x_cg      = 0;  % (ft)
    vtail.sd_x_cg      = 0;  % (ft)
    engn.sd_x_cg       = 0;  % (ft)
    fsys.sd_x_cg       = 0;  % (ft)
    prop.sd_x_cg       = 0;  % (ft)
    payld.sd_x_cg_EOIR = 0;  % (ft)
    payld.sd_x_cg_SAR  = 0;  % (ft)
    payld.sd_x_cg_LiDAR = 0; % (ft)
    payld.sd_x_cg_ANT  = 0;  % (ft)
    payld.sd_x_cg_WR   = 0;  % (ft)
    payld.sd_x_cg_IMU  = 0;  % (ft)
    sfcl.sd_x_cg_wing  = 0;  % (ft)
    sfcl.sd_x_cg_htail = 0;  % (ft)
    sfcl.sd_x_cg_vtail = 0;  % (ft)
end
