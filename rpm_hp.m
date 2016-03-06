%This file is used to find the RPM at various significant velocities.

clear all; clc;

load_base_UAV
load_variation_parameters
load_enviro_parameters
load_unit_conversion
load_airfoils
calc_random_UAV

engn.HP = 3.7; 

W_TO = 37;
hp_rpm = 8500/3.7;
rho_o=0.0023363397;
v_stall = 45;
v_TO = v_stall*1.2;
v_loiter = 50; 
v_cruise = 90;
v_max = 95;

%%
[rho, t, a] = calc_atmos(1000);
v_drag = [(v_stall*mph2fps) (v_loiter*mph2fps)];
M = v_drag/a;

calc_drag;

hp_stall = DRAG.P_t(1);
hp_loiter = DRAG.P_t(2);
rpm_stall = hp_stall*hp_rpm;
rpm_loiter = hp_loiter*hp_rpm;
Cs_stall = (0.638*v_stall*(rho/rho_o)^(1/5))/((hp_stall^(1/5))*(rpm_stall^(2/5)))
Cs_loiter = (0.638*v_loiter*(rho/rho_o)^(1/5))/((hp_loiter^(1/5))*(rpm_loiter^(2/5)))
D_stall =  (v_stall*88)/(rpm_stall*1.6)

%
[rho, t, a] = calc_atmos(7500);
v_drag = [(v_cruise*mph2fps) (v_max*mph2fps)];
M = v_drag/a;

calc_drag;

hp_cruise = DRAG.P_t(1);
hp_max = DRAG.P_t(2);
rpm_cruise = hp_cruise*hp_rpm;
rpm_max = hp_max*hp_rpm;
Cs_cruise=(0.638*v_cruise*(rho/rho_o)^(1/5))/((hp_cruise^(1/5))*(rpm_cruise^(2/5)))
Cs_max=(0.638*v_max*(rho/rho_o)^(1/5))/((hp_max^(1/5))*(rpm_max^(2/5)))
D_max =  (v_max*88)/(rpm_max*1.3)

%%
[rho, t, a] = calc_atmos(0);
v_drag = [v_TO*mph2fps];
M = v_drag/a;

calc_drag;

hp_TO = DRAG.P_t(1);
rpm_TO = hp_TO*hp_rpm;
Cs_TO=(0.638*v_TO*(rho/rho_o)^(1/5))/((hp_TO^(1/5))*(rpm_TO^(2/5)))