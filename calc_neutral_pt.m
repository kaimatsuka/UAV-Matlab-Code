function x_np = calc_neutral_pt(wing, htail, airfoilw, airfoilh)
% calc_neutral_point.m
%
% DESCRIPTION:
%   Calculates the neutral point of the current UAV
%
% INPUT:
%   wing: wing structure
%   htail: horizontal tail structure
%   airfoilw: wing airfoil structure
%   airfoilh: horizontal tail airfoil structure
%
% OUTPUT:
%   x_np: neutral point relative to nose (ft)
%
% REVISION HISTORY:
%   02/27: File created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dereference Block -------------------------------------------------------
chord = wing.c; % Chord length of the wing (ft)
x_w_LE = wing.x_cg - (0.25*chord);  % Distance from nose to wing LE (ft)
x_w_AC = 0.25*chord;    % Distance from wing LE to AC (ft) [quarter chord pt]
x_t_AC = htail.x_cg-wing.x_cg+x_w_AC;   % Distance from wing LE to tail AC (ft)
S_t    = htail.S;   % Horizontal Tail Area (ft2)
S_w    = wing.S;    % Wing Area (ft2)
a_t    = airfoilh.a_t; % Tail Lift Curve Slope (/rad)
a_w    = airfoilw.a_w; % Wing Lift Curve Slope (/rad)

%TODO: CHANGE THIS TO ACTUAL VALUE
ep_alpha = 0.2;     % Downwash Factor

% CALCULATE NEUTRAL POINT RELATIVE TO WING LEADING EDGE
x_np_LE = (x_w_AC+(x_t_AC*(S_t/S_w)*(a_t/a_w)*(1-ep_alpha)))/...
        (1+((S_t/S_w)*(a_t/a_w)*(1-ep_alpha))); % ft
% CALCULATE NEUTRAL POINT RELATIVE TO NOSE
x_np    = x_np_LE + x_w_LE; % ft
