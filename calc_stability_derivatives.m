function SDERIV = calc_stability_derivatives(wing,htail,vtail,DRAG, TRIM)
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
%           M_q             % pitch moment due to pitch rate
%           DRAG.CDi        % induced drag
%           TRIM.i_t        % incident angle
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

load_enviro_parameters

% Parameters
V_H    = (htail.l_T*htail.S)/(wing.S*wing.c); % tail volume ratio
eps_a  = 0.2; % estimated downwash effects
rho    = rho_avg; % density
v_avg  = (V_max+V_stall)/2; % average velocity
v      = V_cruise; % flight velocity
i_t    = -(Cm_ac*CL_a + CM_a*CL)/(CL_a*Cm_i-Cm_a*CL_i); % trim incident angle
eta_h  = 0.5; % TODO: replace with real value
tau_e  = htail.e;
sigma_b = -0.1; % estimated elevator downwash effects

airfoils = airfoil_to_wing(airfoils,wing.A,wing.e);
index    = 1; %TODO: replace with index of chosen airfoil
a_w      = airfoils(index).CL_alpha_rad; % lift curve slope for wing
alpha_0  = airfoils(index).alpha0; 

airfoils = airfoil_to_wing(airfoils,htail.A,htail.e); % lift curve slope for htail
index    = 1; % TODO: replace with index the chosen airfoil
a_t      = airfoils(index).CL_alpha_rad;

% Lift Stability Derivatives-----------------------------------------------

%primary
SDERIV.CL0     =    -a_t*htail.S/wing.S*TRIM.i_t; % depends on airfoil geometry
SDERIV.CL_a    =    a_w+a_t*htail.S/wing.S; % mae 154s lec 10
SDERIV.CL_adot =    2*V_H*SDERIV.CL_a*eps_a; % MAE 154S lec 11
SDERIV.CL_q    =    2*V_H*SDERIV.CL_a; % MAE 154S lec 11
SDERIV.CL_de   =    -tau_e*a_t*htail.S/wing.S; % mae 154s hw4
SDERIV.CL_i    =    -a_t*htail.S/wing.S; % mae 154s lec 10

% Drag Stability Derivatives-----------------------------------------------

SDERIV.CD0     =    DRAG.CDi;            
SDERIV.CD_a    =    DRAG.C_Dairf; % guessed value              
SDERIV.CD_de   =    -SDERIV.CL_a*(htail.l_T*vtail.S)/(wing.S*wing.c)*...
                    tau_e*eta_h*(vtail.b/htail.c);
    % http://faculty.dwc.edu/sadraey/Elevator%20Design.pdf

% Force from Y-Direction---------------------------------------------------

SDERIV.CY_a    =    0; % assume 0 because of symmetric aircraft
SDERIV.CY_beta =    -SDERIV.CY_a*(1-sigma_b)*(v_avg/v).^2*vtail.S/wing.S; 
    % lateral stability pdf
SDERIV.CY_r    =    2*SDERIV.CY_a*(v_avg/v)^2*vtail.S/wing.S/wing.b;
SDERIV.CY_dr   =    SDERIV.CY_de*(v_avg/v)^2*vtail.S/wing.S;

% Roll Moment--------------------------------------------------------------
             
SDERIV.Cl_p    =   -SDERIV.CL_a/12*(1+3*wing.lam)/(1+wing.lam); % mae 154s hw4
SDERIV.Cl_beta =   0; % dihedral term, roll stability
SDERIV.Cl_r    =    SDERIV.CY_r*((vtail.z_cg-z_cg_total)/wing.b*cos(alpha_0)-...
                    (vtail.x_cg-wing.x_cg_total)/wing.b*sin(alpha_0)); % lateral stability pdf
SDERIV.Cl_da   =   2*L/(rho*v^2*wing.S*wing.b*d_a); % lateral stability pdf
SDERIV.Cl_dr   =   SDERIV.CY_dr*(vtail.z_cg-z_cg_total)/wing.b*cos(alpha_0)-...
                    (vtail.x_cg-x_cg_total)/wing.b*sin(alpha_0); % lateral stability pdf

% Pitch Moment-------------------------------------------------------------

SDERIV.Cm0     =   Cm_ac+V_H*i_t; % mae 154s lec 10             
SDERIV.Cm_a    =   -SDERIV.CL_a*htail.l_T/wing.c; % mae154s lec 10               
SDERIV.Cm_adot =   -2*SDERIV.CL_adot*V_H*htail.l_T/wing.C*eps_a; % mae154s lec 11
SDERIV.Cm_q    =   -SDERIV.CL_q*htail.l_T/wing.c;
SDERIV.Cm_de   =   -SDERIV.CL_de*htail.l_T/wing.c;
    % http://faculty.dwc.edu/sadraey/Elevator%20Design.pdf
SDERIV.Cm_i    =   a_t*V_H; % MAE 154s lec 10

% Yaw Moment---------------------------------------------------------------

SDERIV.Cn_beta =    CY_a*(1-sigma_b)*(v_avg/v)^2*vtail.S*htail.l_T/wing.S/wing.b;
SDERIV.Cn_p    =   -0.110; 
SDERIV.Cn_r    =   -CY_r*htail.l_T/wing.b;              
SDERIV.Cn_da   =   2*N/(d_a*rho*v^2*wing.S*wing.b);             
SDERIV.Cn_dr   =   -CY_dr*htail.l_T/wing.b;