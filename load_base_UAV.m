% load_base_UAV.m
%
% DESCRIPTION:
%   This file generate and define the basic parameters of base UAV. 
%
% REVISION HISTORY:
%   02/10: First file is created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Wing --------------------------------------------------------------------

wing.A = 13;      % aspect ratio      
wing.S = 100/13;  % wing area (ft^2) 
wing.lam = 0.5;   % taper ratio (must be between 0 < 1)
wing.lam_q = 0;   % wing quarter chord sweep
% wing.lam_max = 0.49; % sweep of maximum thicknes line 
wing.h = 1.286;   % dist from head to wing 1/4 chord (ft) TODO: where in quarter chord?
wing.Q = 1;       % wing interference factor 
wing.e = 0.8;     % Oswald's efficiency factor 
%wing.gamma = 2;   % wing dihedral (deg)  GLOBAL HAWK

% airfoil properties
%   TODO: update this once we have airfoil for wing
wing.mtr = 0.12; % maximum thickness ratio
wing.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
wing.S_wet = 2.003*wing.S;  % wet area for wing (ft^2)

% Fuselage ----------------------------------------------------------------

fuse.Q = 1.25;    % fuselage interference factor TODO: cite source
fuse.L = 4.2857;  % fuselage max length 
fuse.W = 1.0833;  % fuselage max width 
fuse.D = 1.0833;  % fuselage max depth

% Horizontal Tail ---------------------------------------------------------

% primary
htail.A = 3.9;    % aspect ratio
htail.S = 2.3077; % area (ft^2)
htail.lam = 0.49; % taper ratio of horizontal tail (btw 0 and 1 inclusive)
htail.lam_q = 0;  % horizontal tail sweep angle
% htail.lam_max = 0.49; % sweep of maximum thicknes line 
htail.h = 3.8571; % dist from head to 1/4 chord of horizontal tail (ft)
htail.Q   = 1.08; % horizontal tail interference factor 

% airfoil properties
%   TODO: update this once we have airfoil for horizontal tail
htail.mtr = 0.24; % maximum thickness ratio
htail.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
htail.S_wet = 2.003*htail.S; % wet area for horizontal tail (ft^2)

% Vertical Tail -----------------------------------------------------------

% primay
vtail.S = 1.1539; % area (ft^2)
vtail.A = 1.95;   % aspect ratio (defined as b^2/S) 
vtail.lam = 0.6;  % taper ratio 
vtail.lam_q = 0;  % quarter chord sweep angle
% vtail.lam_max = 0.49; % sweep of maximum thicknes line 
vtail.Q = 1.08;   % interference factor 
vtail.h = 3.8571; % dist from head to 1/4 chord of vertical tail (ft)

% airfoil properties
%   TODO: update this once we have airfoil for vertical tail
vtail.mtr = 0.96; % maximum thickness ratio
vtail.mtl = 0.3;  % chordwise lcoation of the airfoil max thickness location (range 0.3~0.5, Raymer pg 435)
vtail.S_wet = 2.003*vtail.S;     % wet area for vertical tail (ft^2)


% Create baseUAV structure ------------------------------------------------

baseUAV.wing  = wing;
baseUAV.fuse  = fuse;
baseUAV.htail = htail;
baseUAV.vtail = vtail;

