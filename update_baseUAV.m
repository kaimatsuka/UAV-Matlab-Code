function [ new_baseUAV ] = update_baseUAV( UAVsuccess )
% DESCRITPION:
%   This function sets the new base UAV to the last successful UAV with the
%   lightest weight for optimization calculation 
%
% INPUTS:
%   UAVsuccess   = successful UAV structure that's lightest weight
%
% OUTPUTS:
%   new_baseUAV  = has all the basic geometic parameters of the UAV
%
% REVISION HISTORY:
%   2/27 : Function created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wing --------------------------------------------------------------------

new_baseUAV.wing.A = UAVsuccess.wing.A;      % aspect ratio      
new_baseUAV.wing.S = UAVsuccess.wing.S;  % wing area (ft^2) 
new_baseUAV.wing.lam = UAVsuccess.wing.lam;   % taper ratio (must be between 0 < 1)
new_baseUAV.wing.lam_q = UAVsuccess.wing.lam_q;   % wing quarter chord sweep (deg)
% new_baseUAV.wing.lam_max = UAVsuccess.wing.lam_max; % sweep of maximum thicknes line 
new_baseUAV.wing.h_q = UAVsuccess.wing.h_q;   % dist from head to wing 1/4 chord (ft) TODO: where in quarter chord?
%new_baseUAV.wing.gamma = UAVsuccess.wing.gamma;   % wing dihedral (deg)  GLOBAL HAWK

% airfoil properties
new_baseUAV.airfoilw.ind = UAVsuccess.airfoilw.ind;

% Fuselage ----------------------------------------------------------------

new_baseUAV.fuse.L = UAVsuccess.fuse.L;  % fuselage max length 
new_baseUAV.fuse.W = UAVsuccess.fuse.W;  % fuselage max width 
new_baseUAV.fuse.D = UAVsuccess.fuse.D;  % fuselage max depth

% Horizontal Tail ---------------------------------------------------------

% primary
new_baseUAV.htail.A = UAVsuccess.htail.A;    % aspect ratio
new_baseUAV.htail.S = UAVsuccess.htail.S; % area (ft^2)
new_baseUAV.htail.lam = UAVsuccess.htail.lam; % taper ratio of horizontal tail (btw 0 and 1 inclusive)
new_baseUAV.htail.lam_q = UAVsuccess.htail.lam_q;  % horizontal tail sweep angle (deg)
% new_baseUAV.htail.lam_max = UAVsuccess.htail.lam_max; % sweep of maximum thicknes line 
new_baseUAV.htail.h = UAVsuccess.htail.h; % dist from head to 1/4 chord of horizontal tail (ft)

% airfoil properties
new_baseUAV.airfoilh.ind = UAVsuccess.airfoilh.ind;

% Vertical Tail -----------------------------------------------------------

% primay
new_baseUAV.vtail.S = UAVsuccess.vtail.S; % area (ft^2)
new_baseUAV.vtail.A = UAVsuccess.vtail.A;   % aspect ratio (defined as b^2/S) 
new_baseUAV.vtail.lam = UAVsuccess.vtail.lam;  % taper ratio 
new_baseUAV.vtail.lam_q = UAVsuccess.vtail.lam_q;  % quarter chord sweep angle (deg)
% new_baseUAV.vtail.lam_max = UAVsuccess.vtail.lam_max; % sweep of maximum thicknes line 
new_baseUAV.vtail.h = UAVsuccess.vtail.h; % dist from head to 1/4 chord of vertical tail (ft)

% airfoil properties
new_baseUAV.airfoilv.ind = UAVsuccess.airfoilw.ind;

% Propeller ---------------------------------------------------------------

new_baseUAV.prop.h = UAVsuccess.prop.h; % propeller location
new_baseUAV.prop.D = UAVsuccess.prop.D; % propeller diameter (ft)
new_baseUAV.prop.z_cg = UAVsuccess.prop.z_cg;

% Ailron ------------------------------------------------------------------

new_baseUAV.sfcl.ail.S = UAVsuccess.sfcl.ail.S; % area (ft^2) multiplication factor ranges btwn 0 to 0.051
new_baseUAV.sfcl.ail.c = UAVsuccess.sfcl.ail.c;

% Rudder ------------------------------------------------------------------

new_baseUAV.sfcl.rudd.S = UAVsuccess.sfcl.rudd.S; % area (ft^2) multiplication factor ranges btwn 0.3 to 0.5
new_baseUAV.sfcl.rudd.c = UAVsuccess.sfcl.rudd.c;

% Elevator ----------------------------------------------------------------

new_baseUAV.sfcl.elev.S = UAVsuccess.sfcl.elev.S; % area 
new_baseUAV.sfcl.elev.c = UAVsuccess.sfcl.elev.c;

% Fuel --------------------------------------------------------------------

new_baseUAV.fuel = UAVsuccess.fuel;

% Engine ------------------------------------------------------------------
new_baseUAV.engn.ind = UAVsuccess.engn.ind;

% Propeller ---------------------------------------------------------------
new_baseUAV.prop.D = UAVsuccess.prop.D;
new_baseUAV.prop.pitch = UAVsuccess.prop.pitch;

% Payload -----------------------------------------------------------------
new_baseUAV.payld.x_cg_LiDAR = UAVsuccess.payld.x_cg_LiDAR;
new_baseUAV.payld.x_cg_ANT   = UAVsuccess.payld.x_cg_ANT;
new_baseUAV.payld.x_cg_WR    = UAVsuccess.payld.x_cg_WR;
new_baseUAV.payld.x_cg_IMU   = UAVsuccess.payld.x_cg_IMU;

% Weight ----------------------------------------------------------------
new_baseUAV.weight  = UAVsuccess.weight;
end

