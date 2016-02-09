% Optimization code
%
% DESCRIPTION:
%
% REVISION HISTORY:
%   02/09: This file is created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear all; close all;

%% Load Files
load_unit_conversion
load_enviro_parameters
load_requirements

%% Randomize aircraft
uav_params

%% Initial guess of take-off weight 
% W_TO = 45;
W_TO = 30; % initial weight guess of a/c (lbs)
W_tolerance = 0.005; % tolerance of total weight estimate
max_weight_refine = 100; % number of iteration 

%% Refine weight
for ii = 1:max_weight_refine
    
    calc_weight_estimate
    
    if abs(WEIGHT.total-W_TO) < W_tolerance
        break
    else
        W_TO = WEIGHT.total;
    end
    
end

% Check if converged
if ii == max_weight_refine && (abs(WEIGHT.total-W_TO) > W_tolerance)
        error('Weight did not converge')
else
    display(['Weight Converged: W = ', num2str(WEIGHT.total), ' lb']);
    display(['Converged after ', num2str(ii),' iterations']);
end

%% Calculate performance

% TODO:
%     Calculate lift
%     Calculate drag
%     Calculate engine/prop

%% Check performance

% TODO:
%     Check against mission requirement
%     Check stability (CG Calculation)
%     Check volume requirement
%     Check stability derivative

%% Save results

% save aircraft
% save performance
% toss bad UAV, keep good UAV