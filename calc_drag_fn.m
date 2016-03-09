function DRAG = calc_drag_fn(v_drag, alt, W, wing, airfoilw, airfoilh, fuse, htail, vtail, xcg_ttl, static_marg)
%%% calc_drag_fn.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRITPION:
%   This function computes total drag and power required for given
%   velocity. If v_drag is single value, this computes drag for that
%   velocity, and if v_drag is a vector of velocities, this computes the
%   drag for corresponding drag for each velocity.
%
% INPUTS:
%   alt     = altitude (1x1)
%   v_drag  = velcoity of interest(1xM)
%   W       = weight of UAV at this time (1x1)
%   wing    = wing structure
%   airfoilw= airfoil of wing structure
%   fuse    = fuselage structure
%   htail   = horizontal tail structure
%   vtail   = vertical tail structure
%
% OUTPUTS:
%   DRAG.C_L     = lift coefficient (1xM)
%       .C_Lw    = wing component of lift coefficient (1xM)
%       .C_Lh    = tail component of lift coefficient (1xM)
%       .C_Dp    = parasite drag coefficient (1xM)
%       .C_Di    = induced drag coefficient (1xM)
%       .C_Dairf = drag coefficient due to airfoil (1xM)
%       .C_Dt    = total drag coefficient (1xM)
%       .D_p     = parasite drag in lb(1xM)
%       .D_i     = induced drag in lb(1xM)
%       .D_airf  = drag due to airfoil in lb (1xM)
%       .D_h     = drag due to tail in lb (1xM)
%       .D_w     = drag due to wing in lb (1xM)
%       .D_t     = total drag in lb (1xM)
%       .i_t     = incident angle on the tail (1xM)
%       .v       = velocity in ft/s (1xM)
%       .P_p     = power req'd for parasite drag in HP (1xM)
%       .P_i     = power req'd for induced drag in HP (1xM)
%       .P_airf  = power red'd for airfoil drag in HP (1xM)
%       .P_t = total power req'd in HP (1xM)
%
% REVISION HISTORY:
%   02/25: modified from calc_drag.m. Now this is in function form.
%   01/29: fixed form factors equations to use sweep angle of respective
%          components (instead of using that of wing)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_enviro_parameters
load_unit_conversion

% find air density and mach number at particular altitude 
[rho, t, a] = calc_atmos(alt);
M = v_drag/a;

%Parameters (Needed to Calculate Trim Conditions)-------------------------

i_w      = 0; % incident angle of the wing assumed to be 0deg (level flight)
eps_a    = 0.2; % estimated value for downwash effect

alpha_0  = (airfoilw.alpha0)*(pi/180); 
CL       = W./(0.5*rho*v_drag.^2*wing.S); % lift coefficient
alpha    = interp1(airfoilw.CL,airfoilw.alpha,CL)*(pi/180); 
           % finds alpha corresponding to CL value
Cm_ac    = interp1(airfoilw.CL,airfoilw.Cm_ac,CL);
a_w      = airfoilw.a_w;
a_t      = airfoilh.a_t;

l_t = htail.x_cg-xcg_ttl;
V_H      = (l_t*htail.S)/(wing.S*wing.c); 
tau_e    = htail.e; % tail efficiency


%STABILITY DERIVATIVES-----------------------------------------------------

%lift
CL_i    = -a_t*htail.S/wing.S;
CL_de   = tau_e*a_t*htail.S/wing.S;

%drag
Cm_a    = a_w*static_marg;
Cm_i    = a_t*V_H;
Cm_de   = -Cm_i*tau_e;
i_t     = -(Cm_ac*a_w + Cm_a*CL)/(a_w*Cm_i-Cm_a*CL_i); % trim incident angle
Cm0     = Cm_ac+(Cm_i*i_t);
de      = -(Cm0*a_w + Cm_a*CL)/(a_w*Cm_de - Cm_a*CL_de); % elevator deflection angle

%REYNOLDS NUMBER-----------------------------------------------------------

wing.Re   = Re(rho,v_drag,wing.c,mu_a); %average Reynolds number for wing
% winglet.Re = Re(rho,v_drag,winglet.c,mu); %average Reynolds number for winglet
fuse.Re   = Re(rho,v_drag,fuse.L,mu_a); %average Reynolds number for fuselage 
htail.Re  = Re(rho,v_drag,htail.c,mu_a); %average Reynolds number for htail
vtail.Re  = Re(rho,v_drag,vtail.c,mu_a); %average Reynolds number for vtail

%FROM FACTOR--------------------------------------------------------------- 

% TODO:
%   use sweep angle of max thickness line instead of quarter chord
%   use appropriate value for chordwise location of airfoil max thickness ratio (x/c)
%

wing.K    = K(wing.mtl,wing.mtr,M,wing.lam_q); %form factor for wing
% winglet.K = K(winglet.x,winglet.c,winglet.t,M,winglet.lam_q); %form factor for winglet
fuse.K    = ones(size(wing.K))*(1+60/fuse.L^3-fuse.L/400) ; %form factor for fuselage
htail.K   = K(htail.mtl,htail.mtr,M,htail.lam_q); %form factor for horizontal tail  
vtail.K   = K(vtail.mtl,vtail.mtr,M,vtail.lam_q); %form factor for vertical wing

wing.K    = K(wing.mtl,wing.mtr,M,wing.lam_q); %form factor for wing
% winglet.K = K(winglet.x,winglet.c,winglet.t,M,winglet.lam_q); %form factor for winglet
fuse.K    = ones(size(wing.K))*(1+60/fuse.L^3-fuse.L/400) ; %form factor for fuselage
htail.K   = K(htail.mtl,htail.mtr,M,htail.lam_q); %form factor for horizontal tail  
vtail.K   = K(vtail.mtl,vtail.mtr,M,vtail.lam_q); %form factor for vertical wing

%FRICTION COEFFICIENT------------------------------------------------------

wing.c_f    = C_f(wing.Re,M); %skin friction coefficient for wing
% winglet.c_f = C_f(winglet.Re,M); %skin friction coefficient for winglet
fuse.c_f    = C_f(fuse.Re,M); %skin friction coefficient for fuselage
htail.c_f   = C_f(htail.Re,M); %skin friction coefficient for horizontal tail
vtail.c_f   = C_f(vtail.Re,M); %skin friction coefficient for vertical tail

%INTERFERNCE FACTOR--------------------------------------------------------

%  TODO: cite source
wing.Q  = 1;    
fuse.Q  = 1.25;
htail.Q = 1.08;
vtail.Q = 1.08;   

%% Drag Estimation Equations

%MISC DRAG TERMS
C_Dmisc = 0.001;
C_DLP   = 0.001;

% Component drag parameters
% K_vec     = [wing.K; fuse.K; htail.K; vtail.K; winglet.K]; %form factor vector
% Q_vec     = [wing.Q; fuse.Q; htail.Q; vtail.Q; winglet.Q]; %interference factor vector
% C_f_vec   = [wing.c_f; fuse.c_f; htail.c_f; vtail.c_f; winglet.c_f]; %skin friction vector
% S_wet_vec = [wing.S_wet; fuse.S_wet; htail.S_wet; vtail.S_wet; winglet.S_wet]; %wet area vector

% NO WINGLET
K_vec     = [wing.K; fuse.K; htail.K; vtail.K]; %form factor vector
Q_vec     = [wing.Q; fuse.Q; htail.Q; vtail.Q]; %interference factor vector
C_f_vec   = [wing.c_f; fuse.c_f; htail.c_f; vtail.c_f]; %skin friction vector
S_wet_vec = [wing.S_wet; fuse.S_wet; htail.S_wet; vtail.S_wet]; %wet area vector

% Compute lift coefficients
DRAG.C_L     = CL; % Calculated above
DRAG.C_Lh    = (a_t*(((alpha + i_w)*(1-eps_a))-(i_t-i_w)-alpha_0))*(htail.S/wing.S); % lift coefficient contribution of tail
DRAG.C_Lw    = DRAG.C_L - DRAG.C_Lh; % lift coefficient contribution from wing

% Compute drag coefficients
DRAG.C_Dp    = C_Dpi(K_vec,Q_vec,C_f_vec,S_wet_vec,wing.S,C_Dmisc,C_DLP); %parasite dragcoefficient equation
DRAG.C_Diw   = DRAG.C_Lw.*DRAG.C_Lw.*wing.K; % drag contribution of wing
DRAG.C_Dih   = DRAG.C_Lh.*DRAG.C_Lh.*htail.K*htail.S/wing.S; % drag contribution of tail
DRAG.C_Dairf = interp1(airfoilw.CL, airfoilw.Cd, DRAG.C_L); % airfoil produced drag

DRAG.C_Dt = DRAG.C_Dp + DRAG.C_Diw+DRAG.C_Dih+DRAG.C_Dairf;

% Compute drag
DRAG.D_p    = DRAG.C_Dp*0.5*rho.*v_drag.^2*wing.S;
DRAG.D_iw   = DRAG.C_Diw*0.5*rho.*v_drag.^2*wing.S;
DRAG.D_ih   = DRAG.C_Dih*0.5*rho.*v_drag.^2*wing.S;
DRAG.D_airf = DRAG.C_Dairf*0.5*rho.*v_drag.^2*wing.S;

DRAG.D_t    = DRAG.D_p+DRAG.D_iw+DRAG.D_ih+DRAG.D_airf;

% Save incidence angle and 
DRAG.i_t   = i_t;
DRAG.Cm_ac = Cm_ac;
DRAG.v = v_drag;

%POWER---------------------------------------------------------------------

DRAG.P_p    = DRAG.D_p.*v_drag*lbfts2hp;
DRAG.P_iw   = DRAG.D_iw.*v_drag*lbfts2hp;
DRAG.P_ih   = DRAG.D_ih.*v_drag*lbfts2hp;
DRAG.P_airf = DRAG.D_airf.*v_drag*lbfts2hp;
DRAG.P_t    = DRAG.D_t.*v_drag*lbfts2hp;