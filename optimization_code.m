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

% ----- Initialize Counter Variables -----
NUM_ITERATION = 10000;
NUM_SUCCESS = 0;
NUM_FAIL = 0;
NUM_VOLUMEFAILS = 0;
NUM_CLFAILS = 0;
NUM_ENDUFAILS = 0;
NUM_STATICMARGINFAILS = 0;
NUM_LOADFACTFAILS = 0;
NUM_BASEUAVCHANGE = 0;

% Atmos structure has all atmospheric properties for calculation
atmos(1).altitude = 1000;
[atmos(1).rho, atmos(1).T, atmos(1).a] = calc_atmos(atmos(1).altitude);
atmos(2).altitude = 7500;
[atmos(2).rho, atmos(2).T, atmos(2).a] = calc_atmos(atmos(2).altitude);

% ----- INITIAL WEIGHT ESTIMATE OF BASE UAV -----
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%               MONTE CARLO UAV OPTIMIZATION                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Monte Carlo optimization begins.');

for jj = 1:NUM_ITERATION
    
    if mod(jj,10) == 1
        fprintf('.');
    end
    
    if mod(jj,200) == 0
        fprintf('%d iterations\n', jj)  
    end
    
    % ----- SET NEW BASE UAV (GENETIC ALGORITHM) -----
    % Every 2000 iterations, take the lowest successful aircraft weight and
    % if lightwest Monte Carlo AC lighter than base UAV, set that as the
    % new base UAV and optimize for that aircraft
    if (mod(jj,2000) == 0 && NUM_SUCCESS >= 1)
       [w ind] = min(arrayfun(@(x) x.weight.total, UAVsuccess));
       if (w < baseUAV.weight.total)
           baseUAV = update_baseUAV(UAVsuccess(ind));
           NUM_BASEUAVCHANGE = NUM_BASEUAVCHANGE + 1;
       end
    end
    
    % ----- RANDOMIZE AIRCRAFT -----
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                   CALCULATE PERFORMANCE                     %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO:
    %     Calculate engine/prop  NADIA
    %     Calculate moment of inertia
    %     Calculate trim drag(?)
    %     Calculate stability derivatives
    
    % ----- CALCULATE VOLUME REQUIREMENT -----
    length_total = payld.length_TOTAL + fsys.length + engn.length;  % ft
    wing.x_LE = wing.h - (0.25*wing.c);
    htail.x_LE = htail.h - (0.25*htail.c);
    vtail.x_LE = vtail.h - (0.25*vtail.c);
    
    if 0
        v_vec = linspace(V_stall,V_max,100);
        DRAG_1000 = calc_drag_fn(v_vec,atmos(1).altitude,WEIGHT.total,...
                        wing,airfoilw, airfoilh,fuse,htail,vtail); 
        DRAG_7500 = calc_drag_fn(v_vec,atmos(2).altitude,WEIGHT.total,...
                        wing,airfoilw,airfoilh, fuse,htail,vtail);
                    
        figure(1)
        plot(DRAG_1000.v, DRAG_1000.D_t,'b'); hold on, grid on
        plot(DRAG_7500.v, DRAG_7500.D_t,'r');
        legend('1000 ft','7000 ft');
        return 
    end 
    
    
    % ----- CALCULATE LIFT & DRAG -----
    DRAG_1000 = calc_drag_fn([V_stall V_loiter],atmos(1).altitude,WEIGHT.total,...
                        wing,airfoilw, airfoilh,fuse,htail,vtail);
    DRAG_7500 = calc_drag_fn([V_cruise V_max],atmos(2).altitude,...
        WEIGHT.total,wing,airfoilw,airfoilh, fuse,htail,vtail);
    
    % De-reference Lift from DRAG Structure
    CL_max      = DRAG_1000.C_L(1);  % at 1000 ft
    CL_loiter   = DRAG_1000.C_L(2);  % at 1000 ft
    CL_cruise   = DRAG_7500.C_L(1);  % at 7500 ft
    CL_mxspd    = DRAG_7500.C_L(2);  % at 7500 ft
    CL_VEC      = [CL_max CL_loiter CL_cruise CL_mxspd];
    
    
    % ----- CALCULATE ENDURANCE -----
    % Determine how much fuel  burned for climb, cruise, and loiter
    % CLIMBING 
    P_eng_reqd  = [DRAG_1000.P_t DRAG_7500.P_t]/prop.eta_p;
    P_eng_avail = ones(1,4)*engn.HP*prop.eta_p;
    P_excess = (P_eng_avail - P_eng_reqd)*hp2lbfts;    % lbf ft/s
    energy(1) =  ((0.5*WEIGHT.total*V_launch^2)/g)+(WEIGHT.total*5); %lbf ft
    energy(2) =  ((0.5*WEIGHT.total*V_cruise^2)/g)+(WEIGHT.total*atmos(2).altitude); %lbf ft
    endu_climb = (energy(2) - energy(1))/P_excess(2);   % time it takes to climb
    W_fuel(1) = fuel.cp*P_eng_avail(1)*hp2lbfts*endu_climb;
    W_final(1) = WEIGHT.total - W_fuel(1);
    % CRUISING
    [W_final(2), W_fuel(2)] = endu2W_fuel(W_final(1),prop.eta_p,E_max,...
                            fuel.cp, CL_loiter, DRAG_1000.C_Dt(2), atmos(1).rho, wing.S);
    % LOITERING
    [W_final(3), W_fuel(3)] = endu2W_fuel(W_final(2),prop.eta_p,E_min,...
                            fuel.cp, CL_cruise, DRAG_7500.C_Dt(1), atmos(2).rho, wing.S);
    % TOTAL FUEL BURNED
     W_fuel_total = sum(W_fuel);
     
    % ----- CALCULATE LOAD FACTORS ----- 
    % Lift Constrained Load Factor
    loadfact.maxLC_cruise = 0.5*atmos(2).rho*V_cruise^2*(CL_max/(WEIGHT.total/wing.S));
    loadfact.maxLC_loiter = 0.5*atmos(1).rho*V_loiter^2*(CL_max/(WEIGHT.total/wing.S));
    T_avail = (P_avail/V_stall)*hp2lbfts;
    % Thrust Constrained Load Factor
    loadfact.maxTC = max(CL_VEC./[DRAG_1000.C_Dt DRAG_7500.C_Dt])*(T_avail/WEIGHT.total);

    % ----- CALCULATE STATIC MARGIN -----
    calc_CG;    % Center of Mass of UAV
    % Neutral Point
    stab.x_np = calc_neutral_pt(wing, htail, airfoilw, airfoilh); 
    % Static Margin
    stab.static_margin_full = calc_static_marg(stab.x_cg_full,stab.x_np,wing);
    stab.static_margin_empty = calc_static_marg(stab.x_cg_empty,stab.x_np,wing);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Check performance
    %
    % TODO:
    %     Check other various requirements
    %           - volume requirement
                    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%              CHECK AGAINST REQUIREMENTS                     %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FAIL_FLG = 0;
    %---------- Determine if Geometric Fail      -------------
    if(length_total > fuse.L || wing.c+htail.c > fuse.L || wing.x_LE+wing.c > htail.x_LE || ...
            wing.x_LE+wing.c > vtail.x_LE || htail.x_LE+htail.c > fuse.L || ...
            vtail.x_LE+vtail.c > fuse.L)
        NUM_VOLUMEFAILS = NUM_VOLUMEFAILS + 1;
        FAIL_FLG = FAIL_FLG + 1;
        status(FAIL_FLG) = cellstr('Geometric Fail');
    end
    
    %---------- Determine if airfoil provides sufficient lift -------------
    if(CL_max > airfoilw.CLmax) 
       NUM_CLFAILS = NUM_CLFAILS+1;
       FAIL_FLG = FAIL_FLG + 1;
       status(FAIL_FLG) = cellstr('CLmax Fail');
    end
   
   %---------- Determine if we have enough fuel for mission ---------------
   if (W_fuel_total > WEIGHT.fuel)
        NUM_ENDUFAILS = NUM_ENDUFAILS + 1;
        FAIL_FLG = FAIL_FLG + 1;
        status(FAIL_FLG) = cellstr('Fuel Fail');
   end
    
    %--------- Determine if UAV can pull enough G's -----------------------
%     if ((loadfact.maxLC_cruise < N) || (loadfact.maxLC_loiter < N) || (loadfact.maxTC < N))
%         NUM_LOADFACTFAILS = NUM_LOADFACTFAILS + 1;
%         FAIL_FLG = FAIL_FLG + 1;
%         status(FAIL_FLG) = cellstr('Load Factor Fail');
%     end
    
    %---------- Determine if UAV statically stable w/ full & empty fuel ---
    if ((stab.static_margin_full < 0) || (stab.static_margin_empty < 0))
        NUM_STATICMARGINFAILS = NUM_STATICMARGINFAILS + 1;
        FAIL_FLG = FAIL_FLG + 1;
        status(FAIL_FLG) = cellstr('Static Margin Fail');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if FAIL_FLG > 0
        NUM_FAIL = NUM_FAIL + 1;
        UAVfail_ind(NUM_FAIL) = jj;
    else
        NUM_SUCCESS = NUM_SUCCESS + 1;
        UAVsuccess_ind(NUM_SUCCESS) = jj;
        UAVsuccess(NUM_SUCCESS) = saveUAV(wing, airfoilw, fuse, htail, ...
                                airfoilh, vtail, airfoilv, engn, fsys, fuel, prop,...
                                payld, stab, loadfact, WEIGHT, status);
        status(1) = cellstr('Passed Req Check!');
    end
    
    UAVall(jj) = saveUAV(wing, airfoilw, fuse, htail, airfoilh,...
                    vtail, airfoilv, engn, fsys, fuel, prop, payld,...
                    stab, loadfact, WEIGHT, status);
    %%% Save results
    % TODO: if UAV passes criteria, do following
    %     save performance
    %     toss bad UAV, keep good UAV
    %
    % IF UAV PASSES ALL CRITERIA ABOVE:
    % SAVE UAV PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% INPUT DATA
%   history.wing.S
%   ind_good
%   ind_bad

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

% ITERATION PLOTS
figure(2)
subplot(3,6,[1:2 7:8 13:14]), plot_hist(arrayfun(@(x) x.weight.total, UAVall),UAVsuccess_ind,UAVfail_ind), ylabel('Weight (lbs)');
subplot(3,6,3), plot_hist(arrayfun(@(x) x.wing.S, UAVall),UAVsuccess_ind,UAVfail_ind), ylabel('Wing Area (ft^2)');
subplot(3,6,4), plot_hist(arrayfun(@(x) x.wing.A, UAVall),UAVsuccess_ind,UAVfail_ind), ylabel('Wing Aspect Ratio');
subplot(3,6,5), plot_hist(arrayfun(@(x) x.wing.b, UAVall),UAVsuccess_ind,UAVfail_ind), ylabel('Wing Span (ft)');
subplot(3,6,6), plot_hist(arrayfun(@(x) x.wing.c, UAVall),UAVsuccess_ind,UAVfail_ind), ylabel('Wing Chord (ft)');
subplot(3,6,9), plot_hist(arrayfun(@(x) x.htail.S, UAVall),UAVsuccess_ind,UAVfail_ind),ylabel('Horizontal Tail Area (ft^2)');
subplot(3,6,10), plot_hist(arrayfun(@(x) x.htail.A, UAVall),UAVsuccess_ind,UAVfail_ind),ylabel('Horizontal Aspect Ratio');
subplot(3,6,11), plot_hist(arrayfun(@(x) x.htail.b, UAVall),UAVsuccess_ind,UAVfail_ind),ylabel('Horizontal Tail Span (ft)');
subplot(3,6,12), plot_hist(arrayfun(@(x) x.htail.c, UAVall),UAVsuccess_ind,UAVfail_ind),ylabel('Horizontal Tail Chord (ft)');
subplot(3,6,15:18), plot_hist(arrayfun(@(x) x.fuse.L, UAVall),UAVsuccess_ind,UAVfail_ind), ylabel('Fuselage Length (ft)');

% % Plots the lightest aircraft
% [val idx] = min(arrayfun(@(x) x.weight.total, UAVpass));
% figure(1);
% plot_UAV(UAVpass(idx).wing, UAVpass(idx).htail, UAVpass(idx).vtail, UAVpass(idx).fuse, UAVpass(idx).prop);
% 
% figure(3);
% weight_vec = [UAVpass(idx).weight.wing   UAVpass(idx).weight.fuse  ...
%               UAVpass(idx).weight.htail  UAVpass(idx).weight.vtail ...
%               UAVpass(idx).weight.engn   UAVpass(idx).weight.avion ...  
%               UAVpass(idx).weight.fsys   UAVpass(idx).weight.fuel  ...
%               UAVpass(idx).weight.prop];
% weight_labels = {'wing',  'fuselage', 'horizontal tail', 'vertical tail', ...
%                 'engine', 'payload/avionics', 'fuel system', 'fuel',   'propeller'};
% pie(weight_vec, weight_labels); hold on;

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
    
