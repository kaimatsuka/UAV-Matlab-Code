function TRIM = calc_trim_drag_TEMP(CL, htail, wing)
% calc_trim_drag.TEMP.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:
%   This function calculates the the elevator deflection given CL value.
%
% INPUTS:
%      SDERIV.CL_a
%
%
% OUTPUTS:
%         de  % elevator deflection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_enviro_parameters
load_requirements
load_base_UAV
airfoil_analysis

% Parameter

i_w      = 0; % incident angle of the wing assumed to be 0deg (level flight)
eps_a    = 0.2; % estimated value for downwash effect
alpha    = 0;

airfoils = airfoil_to_wing(airfoils,wing.A,wing.e);
index    = 1; % TODO: replace with index of chosen airfoil
alpha_0  = airfoils(index).alpha0;

V_H      = (htail.l_T*htail.S)/(wing.S*wing.c); 
tau_e    = htail.e; % tail efficiency
airfoils = airfoil_to_wing(airfoils,htail.A,htail.e);
index    = 1; % TODO: replace with index of chosen airfoil
a_t      = airfoils(index).CL_alpha_rad; % lift curve slope for tail


% Stability Derivatives-------------------------------------------------

%lift
CL_a = 0; % airfoil property
CL_i = -a_t*htail.S/wing.S;

%drag
Cm_ac = 0; % airfoil property
CL_de = -tau_e*a_t*htail.S/wing.S;
Cm_ac = 0; % airfoil property
Cm_a = CL_a*htail.l_T/wing.c; 
Cm_de = CL_de*htail.l_T/wing.c;
Cm_i = a_t*V_H;

i_t = -(Cm_ac*CL_a + CM_a*CL)/...
            (CL_a*Cm_i-Cm_a*CL_i); % trim incident angle

Cm0 = Cm_ac+V_H*i_t;
        
% TRIM VALUES--------------------------------------------------------------
              
TRIM.de = (Cm0*CL_a + Cm_a*CL)/...
          (CL_a*Cm_de - Cm_a*CL_de); % elevator deflection angle
 
TRIM.CL_t = CL_a((alpha + i_w)*(1-eps_a)+...
            (i_t-i_w)-alpha_0); % lift coefficient of tail

TRIM.CL_w = CL-CL_t; % lift coefficient of wing

