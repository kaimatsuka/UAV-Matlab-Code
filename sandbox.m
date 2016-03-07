% sandbox
% nothing in this file is permanent

clear all, close all, clc
% 
% load_base_UAV
% load_variation_parameters
% calc_random_UAV

% load_base_UAV
load_propeller
load_requirements

prop.D = 1.62;   % propeller diameter (ft)
engn.rpm = 9000; %

V = [10:1:200];
eta_p = calc_propeller(propeller,V,engn.rpm,prop.D);

figure()
plot(V,eta_p);
hold on, grid on
line([V_stall V_stall],[0 0.8],'color','k');
line([V_max   V_max  ],[0 0.8],'color','k');
