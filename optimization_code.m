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
load_requirements

currentPath = pwd;
addpath(genpath(currentPath));

NUM_ITERATION = 10000;
NUM_CLFAILS = 0;
NUM_SUCCESS = 0;
NUM_BASEUAVCHANGE = 0;

% Atmos structure has all atmospheric properties for calculation
atmos(1).altitude = 1000;
[atmos(1).rho, atmos(1).T, atmos(1).a] = calc_atmos(atmos(1).altitude);
atmos(2).altitude = 7500;
[atmos(2).rho, atmos(2).T, atmos(2).a] = calc_atmos(atmos(2).altitude);

display('Monte Carlo optimization begins.');

for jj = 1:NUM_ITERATION
    
    if mod(jj,10) == 1
        fprintf('.');
    end
    
    if mod(jj,200) == 0
        fprintf('%d iterations\n', jj)  
    end
    

    
    % Set new base UAV (genetic algorithm)
    % Compare previous UAV with base UAV; if weight is lower, use that as
    % base
    if (jj > 1 && NUM_SUCCESS >= 1)
       if UAVpass(NUM_SUCCESS).weight.total < baseUAV.weight.total
           baseUAV = UAVpass(NUM_SUCCESS);
           NUM_BASEUAVCHANGE = NUM_BASEUAVCHANGE + 1;
       end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Randomize aircraft
    calc_random_UAV
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    %     Calculate drag
    %     Calculate engine/prop  NADIA
    %     Calculate moment of inertia
    %     Calculate stability (CG Calculation) EUGENE
    %     Calculate trim drag(?)
    %     Calculate stability derivatives
    %
    
    % Calculate Lift
%    V_profile   = [V_stall:1:V_max]; % Velocity Profile for min to max speed
%    [rho, T, a] = calc_atmos((7500+1000)/2); % Avg Values
%    C_L         = (2*W_TO)./(rho*(V_profile.^2)*wing.S); % CLmax is first value
    CL_max      = (2*W_TO)./(atmos(1).rho*V_stall^2*wing.S);  % at 1000 ft
    CL_loiter   = (2*W_TO)./(atmos(1).rho*V_loiter^2*wing.S); % at 1000 ft
    CL_cruise   = (2*W_TO)./(atmos(2).rho*V_cruise^2*wing.S); % at 7500 ft
    CL_mxspd    = (2*W_TO)./(atmos(2).rho*V_max^2*wing.S);    % at 7500 ft
    
    if(CL_max > airfoilw.CLmax) % Determine if airfoil provides sufficient lift
       NUM_CLFAILS = NUM_CLFAILS+1;
       continue                % If not, continue with next iteration
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     calc_CG
    
    % Trim Drag (high altitude scan, 7500 ft)
    M = V_cruise/atmos(2).a;
    DRAG = calc_drag_fn(V_cruise, atmos(2).altitude, W_TO, wing, airfoilw, fuse, htail, vtail);
    TRIMDRAG1 = DRAG;

    % Trim Drag (low altitude loiter, 1000 ft)
    M = V_loiter/atmos(1).a;
    DRAG = calc_drag_fn(V_loiter, atmos(1).altitude, W_TO, wing, airfoilw, fuse, htail, vtail);
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
   
    %%% Save results
    % TODO: if UAV passes criteria, do following
    %     save performance
    %     toss bad UAV, keep good UAV
    %
    % IF UAV PASSES ALL CRITERIA ABOVE:
    % SAVE UAV PARAMETERS
    NUM_SUCCESS = NUM_SUCCESS + 1;
    UAVpass(NUM_SUCCESS).wing     = wing;
    UAVpass(NUM_SUCCESS).airfoilw = airfoilw;
    UAVpass(NUM_SUCCESS).fuse     = fuse;
    UAVpass(NUM_SUCCESS).htail    = htail;
    UAVpass(NUM_SUCCESS).airfoilh = airfoilh;
    UAVpass(NUM_SUCCESS).vtail    = vtail;
    UAVpass(NUM_SUCCESS).airfoilv = airfoilv;
    UAVpass(NUM_SUCCESS).fuse     = fuse;
    UAVpass(NUM_SUCCESS).engn     = engn;
    UAVpass(NUM_SUCCESS).fsys     = fsys;
    UAVpass(NUM_SUCCESS).prop     = prop;
    UAVpass(NUM_SUCCESS).payld    = payld;
    UAVpass(NUM_SUCCESS).weight   = WEIGHT;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     figure(1)
%     plot_UAV(wing,htail,vtail,fuse,prop);

    
end

% INPUT DATA
%   history.wing.S
%   ind_good
%   ind_bad


[sorted_weight sort_ind] = sort(arrayfun(@(x) x.weight.total, UAVpass));; 
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
subplot(3,3,1), plot_hist(arrayfun(@(x) x.weight.total, UAVpass),ind_good,ind_bad), ylabel('Weight (lbs)');
subplot(3,3,2), plot_hist(arrayfun(@(x) x.wing.S, UAVpass),ind_good,ind_bad), ylabel('Wing Area (ft^2)');
subplot(3,3,3), plot_hist(arrayfun(@(x) x.wing.A, UAVpass),ind_good,ind_bad), ylabel('Wing Aspect Ratio');

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
    
