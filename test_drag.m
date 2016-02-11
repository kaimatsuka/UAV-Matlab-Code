% test_drag.m
%
%   Stand alone test code that will plot drag, power required, and power 
%   available as a fucnction of velocity.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all;

% load files
load_unit_conversion
load_requirements
% uav_params
load_UAV_parameters
load_enviro_parameters


%% Initial guess of take-off weight 
W_TO = 46.6; %initial weight guess of a/c (lbs)

%% Input Parameter
v_drag = [60:10:150];      % specify velocity range for plotting
S_ref = wing.S;

% M = 0.0595760168; %Mach number
M = v_drag/a_s; %Mach number

rho_vec = [rho_low rho_avg rho_high];

for ii = 1:length(rho_vec)
    rho = rho_vec(ii);
    calc_drag;
    plot_drag(ii) = DRAG;
end

%% Propeller 

% INPUTS: 
engn.HP = 5.2;

calc_propeller
P_avail = engn.HP*prop.eta_p*ones(1,length(v_drag)); % Power available

%% Plot 

color_set = [0 0   1;  % blue
             1 0   0;  % red
             0 0.5 0]; % green

figure('Name','Drag vs. velocity at 7500 ft')
plot(v_drag,plot_drag(end).D_p(end,:), 'b');  hold on; grid on;
plot(v_drag,plot_drag(end).D_i(end,:), 'r');
plot(v_drag,plot_drag(end).D_t(end,:), 'color', [0 0.5 0]);
line([V_stall V_stall],[0 max(plot_drag(end).D_t(end,:))],'color','k');
line([V_max   V_max  ],[0 max(plot_drag(end).D_t(end,:))],'color','k');
xlabel('Velocity(ft/s)'),ylabel('Drag(lb)');
legend('Parasite','Induced','Total Drag','V_{stall}','V_{max}','location','best');
title('Drag vs. Velcoity');


figure('Name','Total drag vs. velocity at different altitude')
for ii = 1:length(rho_vec)
    plot(v_drag,plot_drag(ii).D_t, 'color', color_set(ii,:));
    hold on; grid on;
end
line([V_stall V_stall],[0 max(plot_drag(end).D_t(end,:))],'color','k');
line([V_max   V_max  ],[0 max(plot_drag(end).D_t(end,:))],'color','k');
xlabel('Velocity(ft/s)'),ylabel('Drag(lb)');
legend('@1000ft','@4150ft','@7500ft','V_{stall}','V_{max}','location','best');
title('Total drag vs. velocity at different altitude');

figure('Name','Power required vs. velocity at 7500 ft')
plot(v_drag,plot_drag(end).P_p(end,:), 'b'); hold on; grid on;
plot(v_drag,plot_drag(end).P_i(end,:), 'r'); 
plot(v_drag,plot_drag(end).P_t(end,:), 'color', [0 0.5 0]); 
line([V_stall V_stall],[0 max(plot_drag(end).P_t(end,:))],'color','k');
line([V_max   V_max  ],[0 max(plot_drag(end).P_t(end,:))],'color','k');
xlabel('Velocity(ft/s)'),ylabel('Power(HP)');
legend('Parasite','Induced','Total Power Req');
title('Power Required vs. Velcoity');


figure('Name','Power required vs. velocity at different altitude')
for ii = 1:length(rho_vec)
    plot(v_drag,plot_drag(ii).P_t, 'color', color_set(ii,:));
    hold on; grid on;
end
plot(v_drag,P_avail, 'm') % power available
line([V_stall V_stall],[0 P_avail(1)],'color','k');
line([V_max   V_max  ],[0 P_avail(1)],'color','k');
xlabel('Velocity(ft/s)'),ylabel('Power(HP)');
legend('@1000ft','@4150ft','@7500ft', 'power available');
title('Power required vs. velocity at different altitude');


figure('Name','Drag Polar')
plot(plot_drag(2).C_Dt,plot_drag(2).C_L); hold on; grid on;
plot(plot_drag(2).C_Dt,afoil.CL_max*ones(size(plot_drag(2).C_Dt)),'r');
xlabel('C_{D}'),ylabel('C_{L}')
title('Drag Polar')
legend('drag polar','max CL')

figure('Name','Cl_p vs. V')
for ii = 1:length(rho_vec)
    plot(v_drag,plot_drag(ii).C_Dp, 'color', color_set(ii,:));
    hold on; grid on;
end
xlabel('Velocity(ft/s)'),ylabel('Parasite Drag Coefficient');
legend('@1000ft','@4150ft','@7500ft');
title('Parasite Drag Coefficient vs. velocity at different altitude');