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

% set path to all subfolders
currentPath = pwd;
addpath(genpath(currentPath));

% ----- Load Files -----
load_unit_conversion
load_enviro_parameters
load_requirements
load_airfoils
load_engine_directory
load_base_UAV
load_variation_parameters


% ----- Initialize Counter Variables -----
NUM_ITERATION = 1000;
NUM_SUCCESS = 0;
NUM_FAIL = 0;
NUM_VOLUMEFAILS = 0;
NUM_CLFAILS = 0;
NUM_ENDUFAILS = 0;
NUM_STATICMARGINFAILS = 0;
NUM_SDERIVFAILS = 0;
NUM_LOADFACTFAILS = 0;
NUM_BASEUAVCHANGE = 0;
NUM_AILFAILS = 0;
NUM_RUDDFAILS = 0;
NUM_ELEVFAILS = 0;


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
        fprintf(' %d iterations\n', jj)  
    end
    
    % ----- SET NEW BASE UAV (GENETIC ALGORITHM) -----
    % Every 2000 iterations, take the lowest successful aircraft weight and
    % if lightwest Monte Carlo AC lighter than base UAV, set that as the
    % new base UAV and optimize for that aircraft
    if (mod(jj,2000) == 0 && NUM_SUCCESS >= 1)
       [w ind] = min(arrayfun(@(x) x.weight.total, UAVsuccess));
       if (w < baseUAV.weight.total)
           baseUAV = update_baseUAV(UAVsuccess(ind));
           calc_non_tunable_parameters;
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
    %     Calculate moment of inertia
    
    % ----- CALCULATE VOLUME REQUIREMENT -----
    length_total = payld.length_TOTAL + fsys.length + engn.length;  % ft
    wing.x_LE = wing.h_q - (0.25*wing.c);
    htail.x_LE = htail.h - (0.25*htail.c);
    vtail.x_LE = vtail.h - (0.25*vtail.c);
    
    % ----- CALCULATE STATIC MARGIN -----
    calc_CG;    % Center of Mass of UAV
    % Neutral Point
    stab.x_np = calc_neutral_pt(wing, htail, airfoilw, airfoilh); 
    % Static Margin
    stab.static_margin_full = calc_static_marg(stab.x_cg_full,stab.x_np,wing);
    stab.static_margin_empty = calc_static_marg(stab.x_cg_empty,stab.x_np,wing);
    
    % ----- CALCULATE PROPELLER EFFICIENCY -----
    calc_propeller
    
    % ----- CALCULATE LIFT & DRAG -----
    if 0
        v_vec = linspace(V_stall,V_max,100);
        DRAG_1000_empty = calc_drag_fn(v_vec,atmos(1).altitude,WEIGHT.total-WEIGHT.fuel,...
                        wing,airfoilw, airfoilh,fuse,htail,vtail,stab.x_cg_empty,stab.static_margin_empty);
        DRAG_1000_full = calc_drag_fn(v_vec,atmos(1).altitude,WEIGHT.total,...
                        wing,airfoilw, airfoilh,fuse,htail,vtail,stab.x_cg_full,stab.static_margin_full);   
        DRAG_7500_empty = calc_drag_fn(v_vec,atmos(2).altitude,WEIGHT.total-WEIGHT.fuel,...
                        wing,airfoilw,airfoilh, fuse,htail,vtail,stab.x_cg_empty,stab.static_margin_empty);
        DRAG_7500_full = calc_drag_fn(v_vec,atmos(2).altitude,WEIGHT.total,...
                        wing,airfoilw,airfoilh, fuse,htail,vtail,stab.x_cg_full,stab.static_margin_full);
       
                    
        figure(1)
        plot(DRAG_1000_empty.v, DRAG_1000_empty.D_t,'b'); hold on, grid on
        plot(DRAG_1000_full.v, DRAG_1000_full.D_t,'-b');
        plot(DRAG_7500_empty.v, DRAG_7500_empty.D_t,'r');
        plot(DRAG_7500_full.v, DRAG_7500_full.D_t,'-r');
        legend('1000 ft Empty','1000 ft Full', '7500 ft Empty', '7500 ft Full');
        return 
    end 
    
    DRAG_1000_empty = calc_drag_fn([V_stall V_loiter],atmos(1).altitude,WEIGHT.total-WEIGHT.fuel,...
                        wing,airfoilw, airfoilh,fuse,htail,vtail,stab.x_cg_empty,stab.static_margin_empty);
    DRAG_1000_full = calc_drag_fn([V_stall V_loiter],atmos(1).altitude,WEIGHT.total,...
                        wing,airfoilw, airfoilh,fuse,htail,vtail,stab.x_cg_full,stab.static_margin_full);
    DRAG_7500_empty = calc_drag_fn([V_cruise V_max],atmos(2).altitude,...
        WEIGHT.total-WEIGHT.fuel,wing,airfoilw,airfoilh, fuse,htail,vtail,stab.x_cg_empty,stab.static_margin_empty);
    DRAG_7500_full = calc_drag_fn([V_cruise V_max],atmos(2).altitude,...
        WEIGHT.total,wing,airfoilw,airfoilh, fuse,htail,vtail,stab.x_cg_full,stab.static_margin_full);
    
    % De-reference Lift from DRAG Structure
    CL_max      = DRAG_1000_full.C_L(1);  % at 1000 ft
    CL_loiter   = DRAG_1000_full.C_L(2);  % at 1000 ft
    CL_cruise   = DRAG_7500_full.C_L(1);  % at 7500 ft
    CL_mxspd    = DRAG_7500_full.C_L(2);  % at 7500 ft
    CL_VEC      = [CL_max CL_loiter CL_cruise CL_mxspd];
    
    
    % ----- CALCULATE ENDURANCE -----
    % Determine how much fuel  burned for climb, cruise, and loiter
    % CLIMBING 
    P_eng_reqd  = [DRAG_1000_full.P_t DRAG_7500_full.P_t]/prop.eta_p;
    P_eng_avail = ones(1,4)*engn.HP*prop.eta_p;
    P_excess = (P_eng_avail - P_eng_reqd)*hp2lbfts;    % lbf ft/s
    energy(1) =  ((0.5*WEIGHT.total*V_launch^2)/g)+(WEIGHT.total*5); %lbf ft
    energy(2) =  ((0.5*WEIGHT.total*V_cruise^2)/g)+(WEIGHT.total*atmos(2).altitude); %lbf ft
    endu_climb = (energy(2) - energy(1))/P_excess(2);   % time it takes to climb
    W_fuel(1) = fuel.cp*P_eng_avail(1)*hp2lbfts*endu_climb;
    W_final(1) = WEIGHT.total - W_fuel(1);
    % CRUISING
    [W_final(2), W_fuel(2)] = endu2W_fuel(W_final(1),prop.eta_p,E_max,...
                            fuel.cp, CL_loiter, DRAG_1000_full.C_Dt(2), atmos(1).rho, wing.S);
    % LOITERING
    [W_final(3), W_fuel(3)] = endu2W_fuel(W_final(2),prop.eta_p,E_min,...
                            fuel.cp, CL_cruise, DRAG_7500_full.C_Dt(1), atmos(2).rho, wing.S);
    % TOTAL FUEL BURNED
     W_fuel_total = sum(W_fuel);
     
    % ----- CALCULATE LOAD FACTORS ----- 
    % Lift Constrained Load Factor
    loadfact.maxLC_cruise = 0.5*atmos(2).rho*V_cruise^2*(CL_max/(WEIGHT.total/wing.S));
    loadfact.maxLC_loiter = 0.5*atmos(1).rho*V_loiter^2*(CL_max/(WEIGHT.total/wing.S));
    T_avail = (engn.HP/V_stall)*hp2lbfts;
    % Thrust Constrained Load Factor
    loadfact.maxTC = max(CL_VEC./[DRAG_1000_full.C_Dt DRAG_7500_full.C_Dt])*(T_avail/WEIGHT.total);
    
    
    % ----- CALCULATE STABILITY DERIVATIVES -----
    SDERIV.alt1000_full = calc_stability_derivatives(atmos(1).rho,V_stall,V_max,V_loiter,wing,airfoilw,...
        htail,airfoilh,vtail,airfoilv,DRAG_1000_full,stab.x_cg_full,stab.z_cg_full,stab.static_margin_full);
    SDERIV.alt1000_empty = calc_stability_derivatives(atmos(1).rho,V_stall,V_max,V_loiter,wing,airfoilw,...
        htail,airfoilh,vtail,airfoilv,DRAG_1000_full,stab.x_cg_empty,stab.z_cg_empty, stab.static_margin_empty);
    SDERIV.alt7500_full = calc_stability_derivatives(atmos(2).rho,V_stall,V_max,V_cruise,wing,airfoilw,...
        htail,airfoilh,vtail,airfoilv,DRAG_7500_full,stab.x_cg_full,stab.z_cg_full, stab.static_margin_full);
    SDERIV.alt7500_empty = calc_stability_derivatives(atmos(2).rho,V_stall,V_max,V_loiter,wing,airfoilw,...
        htail,airfoilh,vtail,airfoilv,DRAG_7500_empty,stab.x_cg_empty,stab.z_cg_empty,stab.static_margin_empty);
    
    % ----- CALCULATE SURFACE CONTROL DEFLECTION -----
    d_a = DRAG_7500_full.C_L/SDERIV.alt7500_full.CL_a*180/pi;
    
    d_r = DRAG_7500_full.C_L/SDERIV.alt7500_full.CL_de*180/pi;
    
    d_e = (-SDERIV.alt7500_full.Cm0*SDERIV.alt7500_full.CL_a+ SDERIV.alt7500_full.Cm_a*DRAG_7500_full.C_L)/...
          (SDERIV.alt7500_full.CL_a*SDERIV.alt7500_full.Cm_de - SDERIV.alt7500_full.Cm_a*SDERIV.alt7500_full.CL_de)*180/pi;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%              CHECK AGAINST REQUIREMENTS                     %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FAIL_FLG = 0;
    clear status
    %---------- Determine if Geometric Fail      -------------
    if(length_total > fuse.L || ...              % PL length > fuselage length
            wing.c+htail.c > fuse.L || ...       % chord lengths of wing & tail > fuselage length
            wing.x_LE+wing.c > htail.x_LE || ... % back of wing is farther back than front of htail
            wing.x_LE+wing.c > vtail.x_LE || ... % back of wing is farther back than front of vtail
            htail.x_LE+htail.c > fuse.L || ...   % back of htail not eq fuselage length
            vtail.x_LE+vtail.c > fuse.L)         % back of vtail > fuselage length
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
    
   %---------- Determine if aileron deflection is maximized ---------------
   if (d_a >= sfcl.ail.maxallow & d_a <= sfcl.ail.minallow)
       NUM_AILFAILS = NUM_AILFAILS + 1;
       FAIL_FLG = FAIL_FLG + 1;
       status(FAIL_FLG) = cellstr('AILERON DEFLECTION EXCEEDED');
   end
   
   %---------- Determine if rudder deflection is maximized ----------------
   
   if (d_r >= sfcl.rudd.maxallow & d_r <= sfcl.rudd.minallow)
      NUM_RUDDFAILS = NUM_RUDDFAILS + 1;
       FAIL_FLG = FAIL_FLG + 1;
       status(FAIL_FLG) = cellstr('RUDDER DEFLECTION EXCEEDED');
   end 
   
   %---------- Determine if elevator deflection is maximized --------------
   if (d_e >= sfcl.elev.maxallow & d_e <= sfcl.elev.minallow)
      NUM_ELEVFAILS = NUM_ELEVFAILS + 1;
       FAIL_FLG = FAIL_FLG + 1;
       status(FAIL_FLG) = cellstr('ELEVATOR DEFLECTION EXCEEDED');
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
    
    %---------- Determine if UAV stability derivatives correct sign -------
    if (check_sderivs(SDERIV.alt1000_full)  == 0 || ...
        check_sderivs(SDERIV.alt1000_empty) == 0 || ...
        check_sderivs(SDERIV.alt7500_full)  == 0 || ...
        check_sderivs(SDERIV.alt7500_empty) == 0)
            NUM_SDERIVFAILS = NUM_SDERIVFAILS + 1;
            FAIL_FLG = FAIL_FLG + 1;
            status(FAIL_FLG) = cellstr('Stability Derivative Fail');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if FAIL_FLG == 0
        NUM_SUCCESS = NUM_SUCCESS + 1;
        UAVsuccess_ind(NUM_SUCCESS) = jj;
        status(1) = cellstr('Passed Req Check!');
        UAVsuccess(NUM_SUCCESS) = saveUAV(wing, airfoilw, fuse, htail, ...
                                airfoilh, vtail, airfoilv, engn, fsys, fuel, prop,...
                                payld, stab, SDERIV, loadfact, sfcl, WEIGHT, status,jj);
    else
        NUM_FAIL = NUM_FAIL + 1;
        UAVfail_ind(NUM_FAIL) = jj;
    end
    
    UAVall(jj) = saveUAV(wing, airfoilw, fuse, htail, airfoilh,...
                    vtail, airfoilv, engn, fsys, fuel, prop, payld,...
                    stab, SDERIV, loadfact, sfcl, WEIGHT, status,jj);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


% if 1
%     figure(10)
%     plot_UAV(wing,htail,vtail,fuse,prop)
% end

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


% Display Monte Carlo Results in Command Prompt
fprintf('\n-----------------------------------------------\n');
disp(['Number of Successful Aircraft: ' num2str(NUM_SUCCESS) ' out of ' num2str(NUM_ITERATION)]);
disp(['  Number of C_Lmax Fails: ' num2str(NUM_CLFAILS)]);
disp(['  Number of Endurance Fails: ' num2str(NUM_ENDUFAILS)]);
disp(['  Number of Load Factor Fails: ' num2str(NUM_LOADFACTFAILS)]);
disp(['  Number of Volume Fails: ' num2str(NUM_VOLUMEFAILS)]);
disp(['  Number of Static Margin Fails: ' num2str(NUM_STATICMARGINFAILS)]);
disp(['  Number of Stability Deriv Fails: ' num2str(NUM_SDERIVFAILS)]);
disp(['  Number of Aileron Deflecion Fails: ' num2str(NUM_AILFAILS)]);
disp(['  Number of Rudder Deflection Fails: ' num2str(NUM_RUDDFAILS)]);
disp(['  Number of Elevator Deflection Fails: ' num2str(NUM_ELEVFAILS)]);


if(NUM_SUCCESS > 0)
    [sorted_weight, sort_ind] = sort(arrayfun(@(x) x.weight.total, UAVsuccess)); 
    ind_light = sort_ind(1:ceil(NUM_SUCCESS/10));     % select good aircraft (lowest 10% in weight)
    ind_heavy  = sort_ind(ceil(NUM_SUCCESS/10)+1:end); % rest of aircraft
    
%     UAV_light_ind = ind_light;
%     UAV_heavy_ind = ind_heavy;
   UAV_light_ind = arrayfun(@(x) x.ind, UAVsuccess(ind_light));
   UAV_heavy_ind = arrayfun(@(x) x.ind, UAVsuccess(ind_heavy));

    [val idx] = min(arrayfun(@(x) x.weight.total, UAVsuccess));
    figure(3);
    plot_UAV(UAVsuccess(idx).wing, UAVsuccess(idx).htail, UAVsuccess(idx).vtail, UAVsuccess(idx).fuse, UAVsuccess(idx).prop);

    figure(4);
    weight_vec = [UAVsuccess(idx).weight.wing   UAVsuccess(idx).weight.fuse  ...
                  UAVsuccess(idx).weight.htail  UAVsuccess(idx).weight.vtail ...
                  UAVsuccess(idx).weight.engn   UAVsuccess(idx).weight.avion ...  
                  UAVsuccess(idx).weight.fsys   UAVsuccess(idx).weight.fuel  ...
                  UAVsuccess(idx).weight.prop];
    weight_labels = {'wing',  'fuselage', 'horizontal tail', 'vertical tail', ...
                    'engine', 'payload/avionics', 'fuel system', 'fuel',   'propeller'};
    pie(weight_vec, weight_labels); hold on;
end


% ITERATION PLOTS -- ALL RUNS
if(NUM_SUCCESS > 0) 
    figure(2)
    % -- WEIGHT HISTOGRAM
    subplot(4,6,[1:2 7:8 13:14 19:20]), 
    plot_hist(arrayfun(@(x) x.weight.total, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind), ylabel('Weight (lbs)');
    % -- WING HISTOGRAMS
    subplot(4,6,3), plot_hist(arrayfun(@(x) x.wing.S, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind), ylabel('Wing Area (ft^2)');
    subplot(4,6,4), plot_hist(arrayfun(@(x) x.wing.A, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind), ylabel('Wing Aspect Ratio');
    subplot(4,6,5), plot_hist(arrayfun(@(x) x.wing.b, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind), ylabel('Wing Span (ft)');
    subplot(4,6,6), plot_hist(arrayfun(@(x) x.wing.c, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind), ylabel('Wing Chord (ft)');
    % -- HORIZONTAL TAIL HISTOGRAMS
    subplot(4,6,9), plot_hist(arrayfun(@(x) x.htail.S, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind),ylabel('Horizontal Tail Area (ft^2)');
    subplot(4,6,10), plot_hist(arrayfun(@(x) x.htail.A, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind),ylabel('Horizontal Aspect Ratio');
    subplot(4,6,11), plot_hist(arrayfun(@(x) x.htail.b, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind),ylabel('Horizontal Tail Span (ft)');
    subplot(4,6,12), plot_hist(arrayfun(@(x) x.htail.c, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind),ylabel('Horizontal Tail Chord (ft)');
    % -- VERTICAL TAIL HISTOGRAMS
    subplot(4,6,15), plot_hist(arrayfun(@(x) x.vtail.S, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind),ylabel('Vertical Tail Area (ft^2)');
    subplot(4,6,16), plot_hist(arrayfun(@(x) x.vtail.A, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind),ylabel('Vertical Aspect Ratio');
    subplot(4,6,17), plot_hist(arrayfun(@(x) x.vtail.b, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind),ylabel('Vertical Tail Span (ft)');
    subplot(4,6,18), plot_hist(arrayfun(@(x) x.vtail.c, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind),ylabel('Vertical Tail Chord (ft)');
    % -- FUSELAGE HISTOGRAMS
    subplot(4,6,21:24), plot_hist(arrayfun(@(x) x.fuse.L, UAVall),UAV_light_ind,UAV_heavy_ind,UAVfail_ind), ylabel('Fuselage Length (ft)');
end


% % ITERATION PLOTS -- ALL RUNS
% if(NUM_SUCCESS > 0) 
%     figure(3)
%     % -- WEIGHT HISTOGRAM
%     subplot(4,6,[1:2 7:8 13:14 19:20]), plot_hist(arrayfun(@(x) x.weight.total, UAVall),UAVsuccess_ind,UAVfail_ind), ylabel('Weight (lbs)');
%     % -- WING HISTOGRAMS
%     subplot(4,6,3), plot_hist(arrayfun(@(x) x.wing.S, UAVall),UAVsuccess_ind,UAVfail_ind), ylabel('Wing Area (ft^2)');
%     subplot(4,6,4), plot_hist(arrayfun(@(x) x.wing.A, UAVall),UAVsuccess_ind,UAVfail_ind), ylabel('Wing Aspect Ratio');
%     subplot(4,6,5), plot_hist(arrayfun(@(x) x.wing.b, UAVall),UAVsuccess_ind,UAVfail_ind), ylabel('Wing Span (ft)');
%     subplot(4,6,6), plot_hist(arrayfun(@(x) x.wing.c, UAVall),UAVsuccess_ind,UAVfail_ind), ylabel('Wing Chord (ft)');
%     % -- HORIZONTAL TAIL HISTOGRAMS
%     subplot(4,6,9), plot_hist(arrayfun(@(x) x.htail.S, UAVall),UAVsuccess_ind,UAVfail_ind),ylabel('Horizontal Tail Area (ft^2)');
%     subplot(4,6,10), plot_hist(arrayfun(@(x) x.htail.A, UAVall),UAVsuccess_ind,UAVfail_ind),ylabel('Horizontal Aspect Ratio');
%     subplot(4,6,11), plot_hist(arrayfun(@(x) x.htail.b, UAVall),UAVsuccess_ind,UAVfail_ind),ylabel('Horizontal Tail Span (ft)');
%     subplot(4,6,12), plot_hist(arrayfun(@(x) x.htail.c, UAVall),UAVsuccess_ind,UAVfail_ind),ylabel('Horizontal Tail Chord (ft)');
%     % -- VERTICAL TAIL HISTOGRAMS
%     subplot(4,6,15), plot_hist(arrayfun(@(x) x.vtail.S, UAVall),UAVsuccess_ind,UAVfail_ind),ylabel('Vertical Tail Area (ft^2)');
%     subplot(4,6,16), plot_hist(arrayfun(@(x) x.vtail.A, UAVall),UAVsuccess_ind,UAVfail_ind),ylabel('Vertical Aspect Ratio');
%     subplot(4,6,17), plot_hist(arrayfun(@(x) x.vtail.b, UAVall),UAVsuccess_ind,UAVfail_ind),ylabel('Vertical Tail Span (ft)');
%     subplot(4,6,18), plot_hist(arrayfun(@(x) x.vtail.c, UAVall),UAVsuccess_ind,UAVfail_ind),ylabel('Vertical Tail Chord (ft)');
%     % -- FUSELAGE HISTOGRAMS
%     subplot(4,6,21:24), plot_hist(arrayfun(@(x) x.fuse.L, UAVall),UAVsuccess_ind,UAVfail_ind), ylabel('Fuselage Length (ft)');
% end
% 
% % Plots the lightest aircraft
% if(NUM_SUCCESS > 0)
%     [sorted_weight, sort_ind] = sort(arrayfun(@(x) x.weight.total, UAVsuccess)); 
%     ind_good = sort_ind(1:ceil(NUM_SUCCESS/10));     % select good aircraft (lowest 10% in weight)
%     ind_heavy  = sort_ind(ceil(NUM_SUCCESS/10)+1:end); % rest of aircraft
%     
%     
%     % ITERATION PLOTS -- SUCCESSFUL RUNS
%     figure(5)
%     % -- WEIGHT HISTOGRAM
%     subplot(4,6,[1:2 7:8 13:14 19:20]), plot_hist(arrayfun(@(x) x.weight.total, UAVsuccess),ind_good,ind_heavy), ylabel('Weight (lbs)');
%     % -- WING HISTOGRAMS
%     subplot(4,6,3), plot_hist(arrayfun(@(x) x.wing.S, UAVsuccess),ind_good,ind_heavy), ylabel('Wing Area (ft^2)');
%     subplot(4,6,4), plot_hist(arrayfun(@(x) x.wing.A, UAVsuccess),ind_good,ind_heavy), ylabel('Wing Aspect Ratio');
%     subplot(4,6,5), plot_hist(arrayfun(@(x) x.wing.b, UAVsuccess),ind_good,ind_heavy), ylabel('Wing Span (ft)');
%     subplot(4,6,6), plot_hist(arrayfun(@(x) x.wing.c, UAVsuccess),ind_good,ind_heavy), ylabel('Wing Chord (ft)');
%     % -- HORIZONTAL TAIL HISTOGRAMS
%     subplot(4,6,9), plot_hist(arrayfun(@(x) x.htail.S, UAVsuccess),ind_good,ind_heavy),ylabel('Horizontal Tail Area (ft^2)');
%     subplot(4,6,10), plot_hist(arrayfun(@(x) x.htail.A, UAVsuccess),ind_good,ind_heavy),ylabel('Horizontal Aspect Ratio');
%     subplot(4,6,11), plot_hist(arrayfun(@(x) x.htail.b, UAVsuccess),ind_good,ind_heavy),ylabel('Horizontal Tail Span (ft)');
%     subplot(4,6,12), plot_hist(arrayfun(@(x) x.htail.c, UAVsuccess),ind_good,ind_heavy),ylabel('Horizontal Tail Chord (ft)');
%     % -- VERTICAL TAIL HISTOGRAMS
%     subplot(4,6,15), plot_hist(arrayfun(@(x) x.vtail.S, UAVsuccess),ind_good,ind_heavy),ylabel('Vertical Tail Area (ft^2)');
%     subplot(4,6,16), plot_hist(arrayfun(@(x) x.vtail.A, UAVsuccess),ind_good,ind_heavy),ylabel('Vertical Aspect Ratio');
%     subplot(4,6,17), plot_hist(arrayfun(@(x) x.vtail.b, UAVsuccess),ind_good,ind_heavy),ylabel('Vertical Tail Span (ft)');
%     subplot(4,6,18), plot_hist(arrayfun(@(x) x.vtail.c, UAVsuccess),ind_good,ind_heavy),ylabel('Vertical Tail Chord (ft)');
%     % -- FUSELAGE HISTOGRAMS
%     subplot(4,6,21:24), plot_hist(arrayfun(@(x) x.fuse.L, UAVsuccess),ind_good, ind_heavy), ylabel('Fuselage Length (ft)');
% 
% end

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