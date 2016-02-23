% Optimization code
%
% DESCRIPTION:
%   This file executes optimization routine of the UAV parameters.
%
% DEPENDENCY:
%   This file calls:
%       load_unit_conversion
%       load_enviro_parameters
%       load_requirements
%       load_base_UAV
%       load_variation_parameters
%       calc_random_UAV
%       calc_weight_estimate
% 
% REVISION HISTORY:
%   02/10: File is populated with pseudo-code. "load_varaition_parameters" 
%          file and "calc_random_UAV" file are added to routine.
%   02/09: This file is created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

%% Load Files
load_unit_conversion
load_enviro_parameters
load_requirements
load_airfoils
load_base_UAV  
load_variation_parameters

currentPath = pwd;
addpath(genpath(currentPath));

NUM_ITERATION = 100;

for jj = 1:NUM_ITERATION
    
    % Set new base UAV (genetic algorithm)
    %
    % TODO: once new good UAV geometries are found, update the base UAV
    %       to new base UAV.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Randomize aircraft
    %
    % TODO: add variability to aircraft parameters, and calculate all
    %       derived geometries.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    calc_random_UAV
    
    % Refine weight

    % W_TO = 45;
    W_TO = 30; % initial weight guess of a/c (lbs)
    W_tolerance = 0.005; % tolerance of total weight estimate
    max_weight_refine = 10; % number of iteration 

    for ii = 1:max_weight_refine

        calc_weight_estimate

        if abs(WEIGHT.total-W_TO) < W_tolerance
            break
        else
            W_TO = WEIGHT.total;
        end

    end

    % Check if converged
    if (ii == max_weight_refine) && (abs(WEIGHT.total-W_TO) > W_tolerance)
        error('Weight did not converge')
    else
        display(['Weight Converged: W = ', num2str(WEIGHT.total), ' lb']);
        display(['Converged after ', num2str(ii),' iterations']);
    end
    
    %%% Calculate performance
    %
    % TODO:
    %     Calculate lift
    %     Calculate drag
    %     Calculate engine/prop  NADIA
    %     Calculate moment of inertia
    %     Calculate stability (CG Calculation) EUGENE
    %     Calculate trim drag(?)
    %     Calculate stability derivatives
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     calc_CG
    
    % Trim Drag (high altitude scan)
    % define inputs needed for calc_drag
    v_drag = 132; %ft/s
    [rho, T, a] = calc_atmos(7500);
    M = v_drag/a;
    W = W_TO;
    calc_drag;
    TRIMDRAG1 = DRAG;

    % Trim Drag (low altitude loiter)
    % define inputs needed for calc_drag
    v_drag = 73.30; %ft/s
    [rho, T, a] = calc_atmos(1000);
    M = v_drag/a;
    W = W_TO;
    calc_drag;
    TRIMDRAG2 = DRAG;
    
    %%% Check performance
    %
    % TODO:
    %     Check against mission requirement
    %     Check other various requirements
    %           - volume requirement
    %           - stability
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     history_wing(jj) = wing;
    
    %%% Save results
    %
    % TODO: if UAV passes criteria, do following
    %     save UAV parameters
    %     save performance
    %     toss bad UAV, keep good UAV
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(1)
    plot_UAV(wing,htail,vtail,fuse,prop);
% 
%     figure(2)
%     hist(history_wing(jj).wing.S)
end

