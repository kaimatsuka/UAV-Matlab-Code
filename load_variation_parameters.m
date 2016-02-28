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

% Wing --------------------------------------------------------------------

wing.sd_A     = 2;  % wing aspect ratio variation std dev 
wing.sd_S     = 2;   % wing area variation (ft^2)
wing.sd_lam   = 0.05; % wing taper ratio variation(n/a)
wing.sd_lam_q = 0;   % wing quarter chord sweep variation (deg)
wing.sd_h     = 0.1;   % head to 1/4 chord distance variation (ft)

% Fuselage ----------------------------------------------------------------

fuse.sd_L     = 0.25; % fuselage max length variation (ft)
fuse.sd_W     = 0; % fuselage max width variation (ft)
fuse.sd_D     = 0; % fuselage max depth variation (ft)

% Horizontal Tail ---------------------------------------------------------

%primary
htail.sd_S       = 0.1; % htail area variation (ft^2)
htail.sd_A       = 0.1; % htail aspect ratio variation std dev
htail.sd_lam     = 0.05; % htail taper ratio variation
htail.sd_lam_q   = 0; % htail quarter chord sweep variation 
% htail.sd_lam_max = 0; % htail sweep of maximum thickness line variation
htail.sd_h       = 0.01; % htail dist from head to 1/4 chord of variation (ft)

% Vertical Tail -----------------------------------------------------------

vtail.sd_S       = 0.1; % vtail area variation (ft^2)
vtail.sd_A       = 0.1; % vtail aspect ratio variation 
vtail.sd_lam     = 0.05; % vtail taper ratio variation
vtail.sd_lam_q   = 0; % vtail quarter chord sweep variation
% vtail.sd_lam_max = 0; % vtail sweep of maximum thickness line variation
vtail.sd_h       = 0.01; % vtail htail interference factor variation

% CG locations ------------------------------------------------------------

wing.sd_x_cg       = 0.5;  % (ft)
fuse.sd_x_cg       = 0.5;  % (ft)
htail.sd_x_cg      = 0.05; % (ft)
vtail.sd_x_cg      = 0.05; % (ft)
engn.sd_x_cg       = 0;  % (ft)
fsys.sd_x_cg       = 0.1;  % (ft)
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