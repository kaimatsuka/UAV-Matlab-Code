% test_drag.m
%
%   Stand alone test code that will plot drag and power required as a 
%   fucnction of velocity.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all;

load_unit_conversion
uav_params
load_enviro_parameters

%% Input Parameter
v_drag = [50:10:200];      % specify velocity range for plotting
S_total = wing.S + htail.S;

% M = 0.0595760168; %Mach number
M = v_drag/a_s; %Mach number

rho_vec = [rho_low rho_avg rho_high];

for ii = 1:length(rho_vec)
    rho = rho_vec(ii);
    calc_drag;
    plot_drag(ii) = DRAG;
end

%% Plot 

color_set = [0 0   1;  % blue
             1 0   0;  % red
             0 0.5 0]; % green

figure('Name','Drag vs. velocity at 7500 ft')
plot(v_drag,plot_drag(end).D_p(end,:), 'b');  hold on; grid on;
plot(v_drag,plot_drag(end).D_i(end,:), 'r');
plot(v_drag,plot_drag(end).D_t(end,:), 'color', [0 0.5 0]);
xlabel('Velocity(ft/s)'),ylabel('Drag(lb)');
legend('Parasite','Induced','Total Drag');
title('Drag vs. Velcoity');


figure('Name','Total drag vs. velocity at different altitude')
for ii = 1:length(rho_vec)
    plot(v_drag,plot_drag(ii).D_t, 'color', color_set(ii,:));
    hold on; grid on;
end
xlabel('Velocity(ft/s)'),ylabel('Drag(lb)');
legend('@1000ft','@4150ft','@7500ft');
title('Total drag vs. velocity at different altitude');

figure('Name','Power required vs. velocity at 7500 ft')
plot(v_drag,plot_drag(end).D_p(end,:), 'b'); hold on; grid on;
plot(v_drag,plot_drag(end).D_i(end,:), 'r'); 
plot(v_drag,plot_drag(end).D_t(end,:), 'color', [0 0.5 0]); 
xlabel('Velocity(ft/s)'),ylabel('Power(HP)');
legend('Parasite','Induced','Total Power Req');
title('Power Required vs. Velcoity');


figure('Name','Power required vs. velocity at different altitude')
for ii = 1:length(rho_vec)
    plot(v_drag,plot_drag(ii).P_t, 'color', color_set(ii,:));
    hold on; grid on;
end
xlabel('Velocity(ft/s)'),ylabel('Power(HP)');
legend('@1000ft','@4150ft','@7500ft');
title('Power required vs. velocity at different altitude');


figure('Name','Drag Polar')
plot(plot_drag(2).C_Dt,plot_drag(2).C_L); hold on; grid on;
xlabel('C_{D}'),ylabel('C_{L}')
title('Drag Polar')