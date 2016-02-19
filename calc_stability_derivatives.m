%% calc_stability_derivatives.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:
%   This function computes the stability derivatives needed to input into
%   the simulation. 
%   
% INPUTS:   
%           l_t             % length from CG to centroid of tail
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
%           lambda          % taper ratio
%           z_v             % z-position of centroid of vtail
%           x_v             % x-position of centroid of vtail
%           z_cg            % z-position of CG
%           x_cg            % x-position of CG
%           v_avg           % average flight velocity
%           v               % flight velocity
%           rho             % density
%           L               % Roll Moment
%           N               % Yaw Moment
%           M_q             % pitch moment due to pitch rate
%           DRAG.CDi        % induced drag
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_H = (l_t*htail.S)/(wing.S*wing.c);

% Lift Stability Derivatives-----------------------------------------------

%primary
SDERIV.CL0     =    0.5; % depends on airfoil geometry
SDERIV.CL_a    =    a_w+a_t*htail.S/wing.S; % mae 154s lec 10
SDERIV.CL_adot =    2*V_H*SDERIV.CL_a*eps_a; % MAE 154S lec 11
SDERIV.CL_q    =    2*V_H*SDERIV.CL_a; % MAE 154S lec 11
SDERIV.CL_de   =    tau_e*a_t*htail.S/wing.S; % mae 154s hw4

% Drag Stability Derivatives-----------------------------------------------

SDERIV.CD0     =    DRAG.CDi;            
SDERIV.CD_a    =    0.430; % guessed value              
SDERIV.CD_de   =    -SDERIV.CL_a*(l_t*vtail.S)/(wing.S*wing.c)*...
                    tau_e*eta_h*(vtail.b/htail.c);
    % http://faculty.dwc.edu/sadraey/Elevator%20Design.pdf

% Force from Y-Direction---------------------------------------------------

SDERIV.CY_a    =    0.5; % guessed value
SDERIV.CY_beta =    -SDERIV.CY_a*(1-sigma_b)*(v_avg/v).^2*vtail.S/wing.S; 
    % lateral stability pdf
SDERIV.CY_r    =    2*SDERIV.CY_a*(v_avg/v)^2*vtail.S/wing.S/wing.b;
SDERIV.CY_dr   =    SDERIV.CY_de*(v_avg/v)^2*vtail.S/wing.S;

% Roll Moment--------------------------------------------------------------

SDERIV.Cl_beta =   -0.023;              
SDERIV.Cl_p    =   -SDERIV.CL_a/12*(1+3*lambda)/(1+lambda); % mae 154s hw4
SDERIV.Cl_r    =    SDERIV.CY_r*((z_v-z_cg)/wing.b*cos(a_0)-...
                    (x_v-x_cg)/wing.b*sin(a_0)); % lateral stability pdf
SDERIV.Cl_da   =   2*L/(rho*v^2*wing.S*wing.b*d_a); % lateral stability pdf
SDERIV.Cl_dr   =   SDERIV.CY_dr*((z_v-z_cg)/wing.b*cos(a_0)-...
                    (x_v-x_cg)/wing.b*sin(a_0)); % lateral stability pdf

% Pitch Moment-------------------------------------------------------------

SDERIV.Cm0     =   Cm_ac+V_H*i_t; % mae 154s lec 10             
SDERIV.Cm_a    =   SDERIV.CL_a*(h_cg-h_n); % mae154s lec 10               
SDERIV.Cm_adot =   -2*SDERIV.CL_adot*V_H*l_t/wing.C*eps_a; % mae154s lec 11
SDERIV.Cm_q    =   1/(rho*v^2*wing.S*wing.c)*M_q; % longitudinal stability pdf
SDERIV.Cm_de   =   -SDERIV.CL_a*eta_h*V_H*vtail.b/htail.b*tau_e;
    % http://faculty.dwc.edu/sadraey/Elevator%20Design.pdf

% Yaw Moment---------------------------------------------------------------

SDERIV.Cn_beta =    CY_a*(1-sigma_b)*(v_avg/v)^2*vtail.S*l_v/wing.S/wing.b;
SDERIV.Cn_p    =   -0.110;
SDERIV.Cn_r    =   -CY_r*l_v/wing.b;              
SDERIV.Cn_da   =   2*N/(d_a*rho*v^2*wing.S*wing.b);             
SDERIV.Cn_dr   =   -CY_dr*l_t/wing.b;
