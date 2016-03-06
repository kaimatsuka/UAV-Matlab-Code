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
new_baseUAV.wing.h = UAVsuccess.wing.h;   % dist from head to wing 1/4 chord (ft) TODO: where in quarter chord?
new_baseUAV.wing.e = UAVsuccess.wing.e;     % Oswald's efficiency factor 
%new_baseUAV.wing.gamma = UAVsuccess.wing.gamma;   % wing dihedral (deg)  GLOBAL HAWK
new_baseUAV.wing.z_cg = UAVsuccess.wing.z_cg;

% airfoil properties
%   TODO: update this once we have airfoil for wing
new_baseUAV.wing.mtr = UAVsuccess.wing.mtr; % maximum thickness ratio
new_baseUAV.wing.mtl = UAVsuccess.wing.mtl;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
new_baseUAV.wing.S_wet = 2.003*UAVsuccess.wing.S;  % wet area for wing (ft^2)

% Fuselage ----------------------------------------------------------------

% fuse.L = 4.2857;  % fuselage max length 
new_baseUAV.fuse.L = UAVsuccess.fuse.L;  % fuselage max length 
new_baseUAV.fuse.W = UAVsuccess.fuse.W;  % fuselage max width 
new_baseUAV.fuse.D = UAVsuccess.fuse.D;  % fuselage max depth
new_baseUAV.fuse.z_cg = UAVsuccess.fuse.z_cg;


% Horizontal Tail ---------------------------------------------------------

% primary
new_baseUAV.htail.A = UAVsuccess.htail.A;    % aspect ratio
new_baseUAV.htail.S = UAVsuccess.htail.S; % area (ft^2)
new_baseUAV.htail.lam = UAVsuccess.htail.lam; % taper ratio of horizontal tail (btw 0 and 1 inclusive)
new_baseUAV.htail.lam_q = UAVsuccess.htail.lam_q;  % horizontal tail sweep angle (deg)
% new_baseUAV.htail.lam_max = UAVsuccess.htail.lam_max; % sweep of maximum thicknes line 
new_baseUAV.htail.h = UAVsuccess.htail.h; % dist from head to 1/4 chord of horizontal tail (ft)
new_baseUAV.htail.z_cg = UAVsuccess.htail.z_cg;

% airfoil properties
%   TODO: update this once we have airfoil for horizontal tail
new_baseUAV.htail.mtr = UAVsuccess.htail.mtr; % maximum thickness ratio
new_baseUAV.htail.mtl = UAVsuccess.htail.mtl;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
new_baseUAV.htail.S_wet = 2.003*UAVsuccess.htail.S; % wet area for horizontal tail (ft^2)

% Vertical Tail -----------------------------------------------------------

% primay
new_baseUAV.vtail.S = UAVsuccess.vtail.S; % area (ft^2)
new_baseUAV.vtail.A = UAVsuccess.vtail.A;   % aspect ratio (defined as b^2/S) 
new_baseUAV.vtail.lam = UAVsuccess.vtail.lam;  % taper ratio 
new_baseUAV.vtail.lam_q = UAVsuccess.vtail.lam_q;  % quarter chord sweep angle (deg)
% new_baseUAV.vtail.lam_max = UAVsuccess.vtail.lam_max; % sweep of maximum thicknes line 
new_baseUAV.vtail.h = UAVsuccess.vtail.h; % dist from head to 1/4 chord of vertical tail (ft)
new_baseUAV.vtail.z_cg = UAVsuccess.vtail.z_cg;


% airfoil properties
%   TODO: update this once we have airfoil for vertical tail
new_baseUAV.vtail.mtr = UAVsuccess.vtail.mtr; % maximum thickness ratio
new_baseUAV.vtail.mtl = UAVsuccess.vtail.mtl;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
new_baseUAV.vtail.S_wet = 2.003*UAVsuccess.vtail.S;     % wet area for vertical tail (ft^2)

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

% Fuel Sys ----------------------------------------------------------------

new_baseUAV.fsys = UAVsuccess.fsys;

% Engine ------------------------------------------------------------------
new_baseUAV.engn.z_cg = UAVsuccess.engn.z_cg;

% Payload -----------------------------------------------------------------
new_baseUAV.payld.z_cg_EOIR  = UAVsuccess.payld.z_cg_EOIR;
new_baseUAV.payld.z_cg_SAR   = UAVsuccess.payld.z_cg_SAR;
new_baseUAV.payld.z_cg_LiDAR = UAVsuccess.payld.z_cg_LiDAR;
new_baseUAV.payld.z_cg_ANT   = UAVsuccess.payld.z_cg_ANT;
new_baseUAV.payld.z_cg_WR    = UAVsuccess.payld.z_cg_WR;
new_baseUAV.payld.z_cg_IMU   = UAVsuccess.payld.z_cg_IMU;

% Weight ----------------------------------------------------------------
new_baseUAV.weight  = UAVsuccess.weight;
end

