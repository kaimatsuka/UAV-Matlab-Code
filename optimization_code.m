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

NUM_ITERATION = 1000;
NUM_SUCCESS = 0;
NUM_CLFAILS = 0;
NUM_ENDUFAILS = 0;
NUM_STATICMARGINFAILS = 0;
NUM_BASEUAVCHANGE = 0;

% Atmos structure has all atmospheric properties for calculation
atmos(1).altitude = 1000;
[atmos(1).rho, atmos(1).T, atmos(1).a] = calc_atmos(atmos(1).altitude);
atmos(2).altitude = 7500;
[atmos(2).rho, atmos(2).T, atmos(2).a] = calc_atmos(atmos(2).altitude);

%INITIAL WEIGHT ESTIMATE OF BASE UAV
% Refine weight
W_TO = 30; % initial weight guess of a/c (lbs)
W_tolerance = 0.005; % tolerance of total weight estimate
max_weight_refine = 10; % number of iteration 

for ii = 1:max_weight_refine

    calc_weight_estimate

    if abs(WEIGHT.total-W_TO) < W_tolerance
        baseUAV.weight = WEIGHT;
        break
    else
        W_TO = WEIGHT.total;
    end

end

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
           baseUAV = update_baseUAV(UAVpass(NUM_SUCCESS));
           NUM_BASEUAVCHANGE = NUM_BASEUAVCHANGE + 1;
       end
    end
    
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
    %     Calculate engine/prop  NADIA
    %     Calculate moment of inertia
    %     Calculate trim drag(?)
    %     Calculate stability derivatives
    %
    
    % Calculate drag
    DRAG_1000 = calc_drag_fn([V_stall V_loiter],atmos(1).altitude,WEIGHT.total,...
                        wing,airfoilw,fuse,htail,vtail);
    DRAG_7500 = calc_drag_fn([V_cruise V_max],atmos(2).altitude,...
        WEIGHT.total,wing,airfoilw,fuse,htail,vtail);
    
    % Calculate Lift  -- Make sure it meets Lift Requirements
    CL_max      = DRAG_1000.C_L(1);  % at 1000 ft
    CL_loiter   = DRAG_1000.C_L(2);  % at 1000 ft
    CL_cruise   = DRAG_7500.C_L(1);  % at 7500 ft
    CL_mxspd    = DRAG_7500.C_L(2);  % at 7500 ft
    
    if(CL_max > airfoilw.CLmax) % Determine if airfoil provides sufficient lift
       NUM_CLFAILS = NUM_CLFAILS+1;
       continue                % If not, continue with next iteration
    end
    
    % Calculate Endurance -- Make sure it meets endurance requirements
    P_eng_reqd  = [DRAG_1000.P_t DRAG_7500.P_t]/prop.eta_p;
    P_eng_avail = ones(1,4)*engn.HP*prop.eta_p;
    P_excess = (P_eng_avail - P_eng_reqd)*hp2lbfts;    % lbf ft/s
    energy(1) =  ((0.5*WEIGHT.total*V_launch^2)/g)+(WEIGHT.total*5); %lbf ft
    energy(2) =  ((0.5*WEIGHT.total*V_cruise^2)/g)+(WEIGHT.total*atmos(2).altitude); %lbf ft
    endu_climb = (energy(2) - energy(1))/P_excess(2);   % time it takes to climb
    W_fuel(1) = fuel.cp*P_eng_avail(1)*hp2lbfts*endu_climb;
    W_final(1) = WEIGHT.total - W_fuel(1);
    [W_final(2), W_fuel(2)] = endu2W_fuel(W_final(1),prop.eta_p,E_max,...
                            fuel.cp, CL_loiter, DRAG_1000.C_Dt(2), atmos(1).rho, wing.S);
    [W_final(3), W_fuel(3)] = endu2W_fuel(W_final(2),prop.eta_p,E_min,...
                            fuel.cp, CL_cruise, DRAG_7500.C_Dt(1), atmos(2).rho, wing.S);
     W_fuel_total = sum(W_fuel);
     if (W_fuel_total > WEIGHT.fuel)    % Determine if we have enough fuel for mission
        NUM_ENDUFAILS = NUM_ENDUFAILS + 1;
        continue;
     end
     
                        
    % Calculate stability (CG + NP Calculation)
    calc_CG;
    stab.x_np = calc_neutral_pt(wing, htail, airfoilw, airfoilh);
    stab.static_margin_full = calc_static_marg(stab.x_cg_full,stab.x_np,wing);
    stab.static_margin_empty = calc_static_marg(stab.x_cg_empty,stab.x_np,wing);
    
    if(stab.static_margin_full < 0) % Determine if UAV statically stable w/ full fuel
        NUM_STATICMARGINFAILS = NUM_STATICMARGINFAILS + 1;
        continue
    elseif(stab.static_margin_empty < 0) % Determine if UAV statically stable w/ empty fuel
        NUM_STATICMARGINFAILS = NUM_STATICMARGINFAILS + 1;
        continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    UAVpass(NUM_SUCCESS).stab     = stab;
    UAVpass(NUM_SUCCESS).weight   = WEIGHT;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     figure(1)
%     plot_UAV(wing,htail,vtail,fuse,prop);

    
end

% INPUT DATA
%   history.wing.S
%   ind_good
%   ind_bad


[sorted_weight sort_ind] = sort(arrayfun(@(x) x.weight.total, UAVpass)); 
ind_good = sort_ind(1:NUM_ITERATION/100);     % select good aircraft (lowest 10% in weight)
ind_bad  = sort_ind(NUM_ITERATION/100+1:end); % rest of aircraft
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
subplot(3,3,4), plot_hist(arrayfun(@(x) x.wing.b, UAVpass),ind_good,ind_bad), ylabel('Wing Span (ft)');
subplot(3,3,5), plot_hist(arrayfun(@(x) x.wing.c, UAVpass),ind_good,ind_bad), ylabel('Wing Chord Length (ft)');
subplot(3,3,6), plot_hist(arrayfun(@(x) x.fuse.L, UAVpass),ind_good,ind_bad), ylabel('Fuselage Length (ft)');


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
    
