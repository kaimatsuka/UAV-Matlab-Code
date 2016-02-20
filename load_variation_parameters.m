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

wing.sd_A     = 0; % wing aspect ratio variation std dev 
wing.sd_S     = 0; % wing area variation (ft^2)
wing.sd_lam   = 0; % wing taper ratio variation
wing.sd_lam_q = 0; % wing quarter chord sweep variation
wing.sd_h     = 0; % head to 1/4 chord distance variation

% Fuselage ----------------------------------------------------------------

fuse.sd_Q     = 0; % fuselage interference factor variation
fuse.sd_L     = 0; % fuselage max length variation
fuse.sd_W     = 0; % fuselage max width variation
fuse.sd_D     = 0; % fuselage max depth variation

% Fuselage ----------------------------------------------------------------

% TODO: populate variation 


% Horizontal Tail ---------------------------------------------------------

% TODO: populate variation 


% Vertical Tail -----------------------------------------------------------

% TODO: populate variation 

% CG locations ------------------------------------------------------------

wing.sd_x_cg = 0;  % (ft)
fuse.sd_x_cg = 0;  % (ft)
htail.sd_x_cg = 0; % (ft)
vtail.sd_x_cg = 0; % (ft)
engn.sd_x_cg = 0;  % (ft)
fsys.sd_x_cg = 0;  % (ft)
prop.sd_x_cg = 0;  % (ft)
payld.sd_x_cg_EOIR = 0;  % (ft)
payld.sd_x_cg_SAR = 0;   % (ft)
payld.sd_x_cg_LiDAR = 0; % (ft)
payld.sd_x_cg_ANT = 0;   % (ft)
payld.sd_x_cg_WR  = 0;   % (ft)
payld.sd_x_cg_IMU = 0;   % (ft)
sfcl.sd_x_cg_wing  = 0;  % (ft)
sfcl.sd_x_cg_htail = 0;  % (ft)
sfcl.sd_x_cg_vtail = 0;  % (ft)