% uav_params.m
%
% DESCRIPTION:
%    Hard coded UAV parameters and their derived parameters
%
% REVISION HISTORY:
%   01/29: moved W_TO to main_script, fix typo for K_i equation. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wing --------------------------------------------------------------------

% primary
wing.A = 13;      %aspect ratio         SCAN EAGLE
wing.S = 100/13;  %wing area (ft^2)     SCAN EAGLE
wing.lam = 0.5;   %taper ratio (must be between 0 < 1)
wing.lam_q = 0;   %wing quarter chord sweep
wing.mtr = 0.12;  %maximum thickness ratio
wing.x = 1.5;     %dist from head to wing 1/4 chord (ft)
wing.Q = 1;       %wing interference factor 
wing.e = 0.8;     %Oswald's efficiency factor 
%wing.gamma = 2;   %wing dihedral (deg)  GLOBAL HAWK

% derived
wing.b = sqrt(wing.A*wing.S);     % wing span (ft)
wing.c = wing.S/wing.b;           % average wing cord (ft)
wing.c_r = 2*wing.c/(1+wing.lam); % wing root chord length (ft)
wing.c_t = wing.c_r*wing.lam;     % wing tip chord length (ft)
wing.t = wing.c*wing.mtr;         % thickness of wing (ft)
wing.S_wet = 2.003*wing.S;        % wet area for wing (ft^2)
wing.K_i = 1/(pi*wing.A*wing.e);  % Drag coefficient


% Winglet -----------------------------------------------------------------

% primary
winglet.Q = 1.08;   %winglet interference factor TODO: where to get this?
winglet.c = 0.001;  %winglet cord length (ft)
winglet.t = 0;      %thickness of winglet (ft)
winglet.S_wet = 0;  %wet winglet area (ft^2) TODO: put some value
winglet.lam_q = 0;  %sweep angle of winglet

% derived
winglet.x = wing.x + 0; %dist from head to winglet 1/4 chord (ft)

% Fuselage ----------------------------------------------------------------

% primary
fuse.Q = 1.25;  %fuselage interference factor TODO: cite source

% derived
fuse.L = wing.b*(3/7);  %fuselage length 
fuse.W = 13/12;         %fuselage max width 
fuse.D = fuse.W;        %fuselage max depth
fuse.r = fuse.L/fuse.W; %fuselage ratio
fuse.S_wet = pi*0.5*((fuse.L*fuse.W)+(fuse.L*fuse.W));%wet fuselage area (ft^2)

% Horizontal Tail ---------------------------------------------------------

% primary
htail.lam = 0.49; %taper ratio of horizontal tail (btw 0 and 1 inclusive)
htail.Q = 1.08; %horizontal tail interference factor 
htail.lam_q = 0;  %horizontal tail sweep angle

% derived
htail.S = 0.3*wing.S; % area (ft)
htail.b = wing.b*0.3; % span (ft) 
htail.t = 0.8*wing.t; % max root thickness (ft)
htail.c = 0.4*wing.c; % average chord length (ft)
htail.x = 0.9*fuse.L; % dist from head to 1/4 chord of horizontal tail (ft)
htail.l_T = htail.x-wing.x; %distance from wing 1/4 MAC to tail 1/4 MAC (ft)
htail.S_wet = 2.003*htail.S; %wet area for horizontal tail 

% Vertical Tail -----------------------------------------------------------

% primay
vtail.c = 5;     % average chord length (ft)
vtail.lam = 0.6; % taper ratio 
vtail.Q = 1.08;  % interference factor 
vtail.lam_q = 0; % quarter chord sweep angle

% derived
vtail.S = htail.S/2;    % area (ft^2)
vtail.b = htail.b/2;    % span (ft)
vtail.t = 0.8*wing.t;   % max root thickness (ft)
vtail.x = 0.9*fuse.L;   % dist from head to 1/4 chord of vertical tail (ft)
vtail.S_wet = 2.003*vtail.S/2; % wet area for vertical tail

% Airfoil -----------------------------------------------------------------

afoil.CL_max = 1.2; % Max Cl