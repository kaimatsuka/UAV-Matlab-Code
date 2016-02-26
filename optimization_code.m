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

clear all; close all; rng(1);

%% Load Files
load_unit_conversion
load_enviro_parameters
load_requirements
load_airfoils
load_base_UAV  
load_variation_parameters

currentPath = pwd;
addpath(genpath(currentPath));

NUM_ITERATION = 10000;

display('Monte Carlo optimization begins.');

for jj = 1:NUM_ITERATION
    
    if mod(jj,10) == 1
        fprintf('.');
    end
    
    if mod(jj,200) == 0
        fprintf('%d iterations\n', jj)  
    end
    
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
    end
    
%         display(['Weight Converged: W = ', num2str(WEIGHT.total), ' lb']);
%         display(['Converged after ', num2str(ii),' iterations']);

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
    DRAG = calc_drag_fn(v_drag, 7500, W, wing, airfoilw, fuse, htail, vtail);
    TRIMDRAG1 = DRAG;

    % Trim Drag (low altitude loiter)
    % define inputs needed for calc_drag
    v_drag = 73.30; %ft/s
    [rho, T, a] = calc_atmos(1000);
    M = v_drag/a;
    W = W_TO;
    DRAG = calc_drag_fn(v_drag, 1000, W, wing, airfoilw, fuse, htail, vtail);
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

    history.wing.A(jj) = wing.A;
    history.wing.S(jj) = wing.S;
    history.weight(jj) = WEIGHT.total;
    
    %%% Save results
    %
    % TODO: if UAV passes criteria, do following
    %     save UAV parameters
    %     save performance
    %     toss bad UAV, keep good UAV
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     figure(1)
%     plot_UAV(wing,htail,vtail,fuse,prop);

    
end

% INPUT DATA
%   history.wing.S
%   ind_good
%   ind_bad

[sorted_weight sort_ind] = sort(history.weight); 
ind_good = sort_ind(1:NUM_ITERATION/10);     % select good aircraft (lowest 10% in weight)
ind_bad  = sort_ind(NUM_ITERATION/10+1:end); % rest of aircraft
% 
% total_data = history.wing.A;
% good_data  = history.wing.A(ind_good);
% bad_data   = history.wing.A(ind_bad);
% 
% max_val = max(total_data);
% min_val = min(total_data);
% order_of_mag = order((max_val-min_val)/10); % order of magnitude
% max_closest = 10^order_of_mag*ceil(max_val/10^order_of_mag);
% min_closest = 10^order_of_mag*floor(min_val/10^order_of_mag);
% edges = [min_closest:10^(order_of_mag+1):max_closest];
% 
% [N_good,BIN] = histc(good_data, edges);
% N_good = fliplr(N_good);
% [N_bad,BIN] = histc(bad_data, edges);
% N_bad = fliplr(N_bad);
% 
% N_total = [N_good; N_bad];
figure(2)
subplot(3,3,1), plot_hist(history.weight,ind_good,ind_bad), ylabel('Weight (lbs)');
subplot(3,3,2), plot_hist(history.wing.S,ind_good,ind_bad), ylabel('Wing Area (ft^2)');
subplot(3,3,3), plot_hist(history.wing.A,ind_good,ind_bad), ylabel('Wing Aspect Ratio');

% figure(2)
% hist(history.wing.A,10), hold on, 
%     subplot(3,3,1)
%     hist(history.wing.A(ind_good))
%     title('Wing aspect ratio')
%     subplot(3,3,2)
%     hist(history.wing.S,10)
%     title('Wing area'), xlabel('ft^{2}')
%     subplot(3,3,3)
%     hist(history.weight,10)
%     title('Total Weight'), xlabel('lb')
%     subplot(3,3,4)
%     [n1, xout1] = hist(history.weight(ind_good));
%     bar(xout1,n1,'r'); grid on, hold on,
%     [n2, xout2] = hist(history.weight(ind_bad));
%     bar(xout2,n2,'g');
    
