% uav_params.m
%
% DESCRIPTION:
%    Hard coded UAV parameters and their derived parameters
%
% REVISION HISTORY:
%   01/29: moved W_TO to main_script, fix typo for K_i equation. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PRIMARY UAV Parameters
%
%  They are hard coded numbers, fundamental parameters of UAV.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e_H = 0.8;

% Wing --------------------------------------------------------------------

% primary
wing.b   = 10;      %wing span(ft)
wing.mtr = 0.12;    %maximum thickness ratio
wing.swp_ang = 0;   %sweep angle %TODO: is this same thing as lam_q?
wing.lam = 2;       %taper ratio
wing.lam_q = 0;     %wing quarter chord sweep
wing.x = 1.5;       %dist from head to wing 1/4 chord (ft)
wing.Q = 1;         %wing interference factor 
wing.e = 0.8;       %Oswald's efficiency factor 

% derived
wing.c = wing.b/13; %average wing chord length(ft)
wing.t = wing.c*wing.mtr; %thickness of wing (ft)
wing.S = wing.b*wing.c; %wing area(ft^2)(reference area)
wing.A = (wing.b^2)/wing.S; %aspect ratio
wing.S_wet = 2.003*wing.S; %wet area for wing (ft^2)
wing.K_i = 1/(pi*wing.A*wing.e);  %Drag coefficient

% Winglet -----------------------------------------------------------------

%primary
winglet.c = 0.001;  %winglet cord length (ft)
winglet.t = 0;      %thickness of winglet (ft)
winglet.x = 1.5;    %dist from head to winglet 1/4 chord (ft)
winglet.Q = 1.08;   %winglet interference factor TODO: where to get this?
winglet.S_wet = 0;  %wet winglet area (ft^2) TODO: put some value

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

htail.S = 0.3*wing.S; %horizontal tail area (ft)
htail.b = wing.b*0.3; %horizontal tail span (ft) 
htail.t = 0.8*wing.t; %horizontal tail max root thickness (ft)
htail.c = 0.4*wing.c; %average horizontal tail chord length (ft)
htail.x = 0.9*fuse.L; %(ft)
htail.l_T = htail.x-wing.x; %distance from wing 1/4 MAC to tail 1/4 MAC (ft)
htail.S_wet = 2.003*htail.S; %wet area for horizontal tail 
htail.lam = 0.49; %taper ratio of horizontal tail
htail.Q = 1.08; %horizontal tail interference factor 

% Vertical Tail -----------------------------------------------------------

vtail.S = htail.S/2;    % vertical tail area (ft^2)
vtail.b = htail.b/2;    % vertical tail span (ft)
vtail.c = 5;            % average vertical tail chord length (ft)
vtail.t = 0.8*wing.t;   % vertical tail max root thickness (ft)
vtail.x = 0.9*fuse.L;
vtail.Q = 1.08;         % vertical tail interference factor 
vtail.S_wet = 2.003*vtail.S/2; % wet area for vertical tail
vtail.lam = 0.6; % taper ratio 
