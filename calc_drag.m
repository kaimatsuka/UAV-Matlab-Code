%%% calc_drag.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRITPION:
%   This function computes total drag and power required for given
%   velocity. If v_drag is single value, this computes drag for that
%   velocity, and if v_drag is a vector of velocities, this computes the
%   drag for corresponding drag for each velocity.
%
% INPUTS:
%   rho     = density of air (1x1)
%   v_drag  = velcoity of interest(1xM)
%   S_total = total lift surface area (1x1)
%   M       = Mach number (1xM)
%
% OUTPUTS:
%   DRAG.C_L  = lift coefficient (1xM)
%       .C_Dp = parasite drag coefficient (1xM)
%       .C_Di = induced drag coefficient (1xM)
%       .C_Dt = total drag coefficient (1xM)
%       .D_p = parasite drag in lb(1xM)
%       .D_i = induced drag in lb(1xM)
%       .D_t = total drag in lb (1xM)
%       .v   = velocity in ft/s (1xM)
%       .P_p = power req'd for parasite drag in HP (1xM)
%       .P_i = power req'd for induced drag in HP (1xM)
%       .P_t = total power req'd in HP (1xM)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%REYNOLDS NUMBER-----------------------------------------------------------

wing.Re   = Re(rho,v_drag,wing.c,mu); %average Reynolds number for wing
winglet.Re = Re(rho,v_drag,winglet.c,mu); %average Reynolds number for winglet
fuse.Re   = Re(rho,v_drag,fuse.L,mu); %average Reynolds number for fuselage 
htail.Re  = Re(rho,v_drag,htail.c,mu); %average Reynolds number for htail
vtail.Re  = Re(rho,v_drag,vtail.c,mu); %average Reynolds number for vtail

%FROM FACTOR--------------------------------------------------------------- 

wing.K    = K(wing.x,wing.c,wing.t,M,wing.swp_ang); %form factor for wing
winglet.K = K(winglet.x,winglet.c,winglet.t,M,wing.swp_ang); %form factor for winglet
fuse.K    = ones(size(wing.K))*(1+60/fuse.L^3-fuse.L/400) ; %form factor for fuselage
htail.K   = K(htail.x,htail.c,htail.t,M,wing.swp_ang); %form factor for horizontal tail  
vtail.K   = K(vtail.x,vtail.c,vtail.t,M,wing.swp_ang); %form factor for vertical wing

%FRICTION COEFFICIENT------------------------------------------------------

wing.c_f    = C_f(wing.Re,M); %skin friction coefficient for wing
winglet.c_f = C_f(winglet.Re,M); %skin friction coefficient for winglet
fuse.c_f    = C_f(fuse.Re,M); %skin friction coefficient for fuselage
htail.c_f   = C_f(htail.Re,M); %skin friction coefficient for horizontal tail
vtail.c_f   = C_f(vtail.Re,M); %skin friction coefficient for vertical tail

%% Drag Estimation Equations

%MISC DRAG TERMS
C_Dmisc = 0.001;
C_DLP   = 0.001;

% Component drag parameters
K_vec     = [wing.K; fuse.K; htail.K; vtail.K; winglet.K]; %form factor vector
Q_vec     = [wing.Q; fuse.Q; htail.Q; vtail.Q; winglet.Q]; %interference factor vector
C_f_vec   = [wing.c_f; fuse.c_f; htail.c_f; vtail.c_f; winglet.c_f]; %skin friction vector
S_wet_vec = [wing.S_wet; fuse.S_wet; htail.S_wet; vtail.S_wet; winglet.S_wet]; %wet area vector


% Compute drag coefficients
DRAG.C_L  = W_TO./(0.5*rho*v_drag.^2*S_total);
DRAG.C_Dp = C_Dpi(K_vec,Q_vec,C_f_vec,S_wet_vec,wing.S,C_Dmisc,C_DLP); %parasite dragcoefficient equation
DRAG.C_Di = wing.K_i*DRAG.C_L.^2; %induced drag coefficient equation
DRAG.C_Dt = DRAG.C_Dp + DRAG.C_Di;

DRAG.D_p = DRAG.C_Dp*0.5*rho.*v_drag.^2*S_total;
DRAG.D_i = DRAG.C_Di*0.5*rho.*v_drag.^2*S_total;
DRAG.D_t = DRAG.D_p+DRAG.D_i;

DRAG.v = v_drag;

%POWER---------------------------------------------------------------------

DRAG.P_p = DRAG.D_p.*v_drag*lbfts2hp;
DRAG.P_i = DRAG.D_i.*v_drag*lbfts2hp;
DRAG.P_t = DRAG.D_t.*v_drag*lbfts2hp;
