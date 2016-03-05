% cacl_weight_estimate.m
%
% DESCRIPTION:
%   This funciton assumes that aircraft parameters are already loaded.
%
% INPUT:
%   W_TO = initial guess of weight at take off
%
% REVISION HISTORY:
%   02/09: File created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Weight Estimation

% Nicoli Weight Estimat Equation
WEIGHT.wing  = W_w(wing.lam, wing.S, W_TO, N, wing.A, wing.lam_q, wing.mtr, V_e); %wing 
WEIGHT.fuse  = W_f(W_TO,N,fuse.L,fuse.W,fuse.D,V_e); %fuselage
WEIGHT.htail = W_ht(W_TO,N,htail.l_T, htail.S, htail.b, htail.t_r); %horizontal tail
WEIGHT.vtail = W_vt(W_TO,N,vtail.S,vtail.b,vtail.t_r); %vertical tail
WEIGHT.fsys  = W_FS(fuel.V,fsys.int,fsys.N_T, engn.N_E); %fuel system %TODO: update with actual fuel tank 
% WEIGHT.engn  = W_p(engn.W_bare,engn.N_E); %propulsion
% WEIGHT.avion = W_TRON(W_AU); %electronics/avionics
% WEIGHT.sc    = W_SC(W_TO); %surface controls (powered)
% WEIGHT.esys  = W_ES(WEIGHT.fsys,WEIGHT.avion); %electrical system

% Sum of components 
WEIGHT.engn  = 1.1*engn.W_bare*engn.N_E; % assume 10% for mounting fudge factor
WEIGHT.avion = payld.w_total;  % electronics/avionics
WEIGHT.fuel  = fuel.W;
WEIGHT.prop  = 0.75; %TODO: incorporate equation

% Detailed weight list
w_detail_vec = [WEIGHT.wing;
                WEIGHT.fuse;
                WEIGHT.htail;
                WEIGHT.vtail;
                WEIGHT.engn;
                WEIGHT.fsys;
                WEIGHT.prop; 
                WEIGHT.fuel;
                payld.w_EOIR;
                payld.w_SAR;
                payld.w_LiDAR;
                payld.w_ANT;
                payld.w_IMU;
                payld.w_WR;
                sfcl.ail.sc_W;
                sfcl.rudd.sc_W;
                sfcl.elev.sc_W];

WEIGHT.total = sum(w_detail_vec);
