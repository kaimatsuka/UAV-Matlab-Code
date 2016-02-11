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


%% Engine sizing calculation

% calc_engn

engn.W_bare = 3.5; % bare engine weight
engn.N_E = 1;      % number of engines
engn.HP = 5.2;     %engine horse power (selected from engine)

%% Propeller sizing calcuation

calc_propeller
P_avail = engn.HP*prop.eta_p; % Power available

%% Fuel Calculation

fuel.cp = 0.85/550/3600;  %[1/ft] specific fuel consumption
fuel.rho = 6.073; %[lbm/gallon] density of fuel for octane gas
fuel.W = 3.5;             %[lb] fuel weight (calculated using test_fuel.m file)
fuel.V = fuel.W/fuel.rho; %[gallon] volume of fuel
fuel.V = fuel.V*gallon2ft3; %[ft^3] volume of fuel

%% Fuel Calculation

fsys.int = 0.00; %percent of fuel tanks that are integral TODO: find out how much is integral
fsys.N_T = 1;    %number of separate fuel tanks

%% Payload 

payld.w_EOIR  = 12.57; % EO/IR weight(lbs)
payld.w_SAR   = 2;     % Synthetic Apperature Radar weight(lbs)
payld.w_LiDAR = 1;     % LiDAR weight(lbs)
payld.w_ANT   = 0.15;  % UHF/VHF antenna weight (lbs)
payld.w_WR    = 0.25;  % WaveRelay weight (lbs)
payld.w_IMU   = 0.11;  % Ellipse E mini (lbs)

% Pyaload/Avionics total weight
payld.w_total = payld.w_EOIR + payld.w_SAR + payld.w_LiDAR + payld.w_ANT ...
                + payld.w_WR + payld.w_IMU;

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
WEIGHT.engn  = 1.1*engn.W_bare*engn.N_E; % 10% mounting fudge factor
WEIGHT.avion = payld.w_total;  % electronics/avionics
WEIGHT.sc    = 0.4915*W_TO^(2/3); % Nadia Equation %TODO: update with actual servos
WEIGHT.fuel  = fuel.W;
WEIGHT.prop  = 0.75; %TODO: incorporate equation

% Detailed weight list
w_detail_vec = [WEIGHT.wing ...
                WEIGHT.fuse ...
                WEIGHT.htail ...
                WEIGHT.vtail ...
                WEIGHT.engn ...
                WEIGHT.fsys ...
                WEIGHT.prop ... 
                WEIGHT.fuel ...
                payld.w_EOIR ...
                payld.w_SAR ...
                payld.w_LiDAR ...
                payld.w_ANT ...
                payld.w_IMU ...
                payld.w_WR ...
                WEIGHT.sc*0.75 ...
                WEIGHT.sc*0.125 ...
                WEIGHT.sc*0.125];

WEIGHT.total = sum(w_detail_vec);
