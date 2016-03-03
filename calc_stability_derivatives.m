function SDERIV = calc_stability_derivatives(rho, V_stall,V_max,V_cruise,wing,airfoilw,...
    htail,airfoilh,vtail,airfoilv,DRAG,x_cg_total,z_cg_total, static_margin)
%% calc_stability_derivatives.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:
%   This function computes the stability derivatives needed to input into
%   the simulation. 
%   
% INPUTS:   
%           htail.l_T       % length from CG to centroid of tail
%           a_t             % lift slope for tail
%           a_w             % lift slope for wing
%           a_0             % zero angle of attack
%           wing.S          % wing area
%           wing.c          % wing chord
%           wing.b          % wing span
%           htail.S         % htail area
%           htail.c         % htail chord
%           htail.b         % htail span
%           vtail.S         % vtail area
%           vtail.c         % vtail chord
%           vtail.b         % vtail span
%           tau_e           % effectiveness of AoA
%           eta_h           % htail dynamic pressure ratio
%           sigma_b         % downwash effect on elevator (negative)
%           eps_a           % downwash effect
%           d_a             % aeliron deflection
%           wing.lam        % taper ratio
%           vtail.z_cg      % z-position of CG of vtail
%           vtail.x_cg      % x-position of CG of vtail
%           z_cg_total      % z-position of CG
%           x_cg_total      % x-position of CG
%           v_avg           % average flight velocity
%           v               % flight velocity
%           rho             % density
%           L               % Roll Moment
%           N               % Yaw Moment
%           DRAG.CDi        % induced drag
%           DRAG.i_t        % incident angle on the tail
%
% OUTPUTS:
%   Lift: 
%    SDERIV.CL0
%           CL_a
%           CL_adot
%           CL_q 
%           CL_de
%   Drag:
%    SDERIV.CD0
%           CD_a
%           CD_de
%   Pitch Moments:
%    SDERIV.Cm0
%           Cm_a
%           Cm_adot
%           Cm_q
%           Cm_de
%   Force in Y-direction:
%    SDERIV.CY_beta
%           CY_dr
%   Roll Moments:
%    SDERIV.Cl_beta
%           Cl_p
%           Cl_r
%           Cl_da
%           Cl_dr
%   Yaw Moments:
%    SDERIV.Cn_beta
%           Cn_p
%           Cn_r
%           Cn_da
%           Cn_dr
% 
% Revision History
%   02/12: function created
%   02/18: input correct equations for the coefficients
%   02/27: change the variables to match code variable names
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_mom = 0.5;    % TODO: CHANGE TO REAL VALUE
d_a = 0.1;      % TODO: CHANGE TO REAL VALUE
Cm_ac = DRAG.Cm_ac;

% Parameters
l_t = htail.x_cg-x_cg_total;
l_v = vtail.x_cg-x_cg_total;
V_H    = (l_t*htail.S)/(wing.S*wing.c); % tail volume ratio
V_V    = (l_v*vtail.S)/(wing.S*wing.b); % vertical tail volume 
eps_a  = 0.2; % estimated downwash effects
v      = V_cruise; % flight velocity
eta_h  = 0.5; % TODO: replace with real value
tau_e  = htail.e;
sigma_b = -0.1; % estimated elevator downwash effects

a_w      = airfoilw.a_w; % lift curve slope for wing
alpha_0  = (airfoilw.alpha0)*(pi/180); 

a_t      = airfoilh.a_t;   % /rad
a_v      = airfoilv.a_v;    % /rad

% Lift Stability Derivatives-----------------------------------------------

%primary
SDERIV.CL0     =    -a_t*htail.S/wing.S*DRAG.i_t; % depends on airfoil geometry
SDERIV.CL_a    =    a_w+(a_t*(htail.S/wing.S)*(1-eps_a)); % mae 154s lec 10
SDERIV.CL_adot =    2*V_H*a_t*eps_a; % MAE 154S lec 11
SDERIV.CL_q    =    2*V_H*a_t; % MAE 154S lec 11
SDERIV.CL_de   =    tau_e*a_t*htail.S/wing.S; % mae 154s hw4
SDERIV.CL_i    =    -a_t*htail.S/wing.S; % mae 154s lec 10

% Drag Stability Derivatives-----------------------------------------------

SDERIV.CD0     =    DRAG.C_Dp;            
SDERIV.CD_a    =    DRAG.C_Dairf; % guessed value              
SDERIV.CD_de   =    -SDERIV.CL_a*((l_t*vtail.S)/(wing.S*wing.c))*...
                    tau_e*eta_h*(vtail.b/htail.c);
    % http://faculty.dwc.edu/sadraey/Elevator%20Design.pdf

% Force from Y-Direction---------------------------------------------------

SDERIV.CY_de   =    0.1; %TODO: CHANGE THIS
SDERIV.CY_beta =    -vtail.S*a_v/wing.S; % McCormick Eqn 10.57
SDERIV.CY_p    =    -vtail.S*a_v*(vtail.z_cg-wing.z_cg); % McCormick Eqn 10.59
SDERIV.CY_r    =    vtail.S*a_v*l_v/wing.S; % McCormick Eqn 10.58
SDERIV.CY_dr   =    vtail.S*a_v*tau_e/wing.S;   % McCormick Eqn 9.127

% Roll Moment--------------------------------------------------------------
             
SDERIV.Cl_p    =   -SDERIV.CL_a/12*(1+3*wing.lam)/(1+wing.lam); % mae 154s hw4
SDERIV.Cl_beta =   0; % dihedral term, roll stability
SDERIV.Cl_r    =   ((SDERIV.CL_a/6)*((1+(3*wing.lam))/(1+wing.lam)));  % McCormick Eqn 9.128
SDERIV.Cl_dr   =   SDERIV.CY_dr*(vtail.z_cg-z_cg_total)/wing.b*cos(alpha_0)-...
                    (vtail.x_cg-x_cg_total)/wing.b*sin(alpha_0); % lateral stability pdf

% Pitch Moment-------------------------------------------------------------

SDERIV.Cm0     =   Cm_ac+(V_H*a_t*DRAG.i_t); % mae 154s lec 10             
SDERIV.Cm_a    =   -SDERIV.CL_a*static_margin; % mae154s lec 10               
SDERIV.Cm_adot =   -2*a_t*V_H*l_t*eps_a/wing.c; % mae154s lec 11
SDERIV.Cm_q    =   -2*a_t*l_t*V_H/wing.c;
SDERIV.Cm_de   =   -SDERIV.CL_de*htail.l_T/wing.c;
    % http://faculty.dwc.edu/sadraey/Elevator%20Design.pdf
SDERIV.Cm_i    =   a_t*V_H; % MAE 154s lec 10

% Yaw Moment---------------------------------------------------------------

SDERIV.Cn_beta =    V_V*a_v*(1-eps_a);  % McCormick Eqn 9.114
SDERIV.Cn_p    =   -0.110; 
%SDERIV.Cn_r    =   -SDERIV.CY_r*l_t/wing.b;              
SDERIV.Cn_da   =   2*N_mom/(d_a*rho*v^2*wing.S*wing.b);             
SDERIV.Cn_dr   =   -SDERIV.CY_dr*l_t/wing.b;