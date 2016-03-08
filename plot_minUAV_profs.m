function [ ] = plot_minUAV_profs(UAVmin, atmos, V_stall, V_loiter, V_cruise, V_max, propeller)
%% plotDragProf.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:
%   This function plots the drag vs velocity and power vs velocity profiles
%   for the lightest weight aircraft
%   
% INPUTS:
%   UAVstruct: structure of lightest optimized UAV
%   atmos : atmospheric properties
%   V_stall : stall speed (ft/s)
%   V_loiter : loiter speed (ft/s)
%   V_cruise : cruise speed (ft/s)
%   V_max   : max speed (ft/s)
% 
% OUTPUTS:
%   figure w/ drag vs. velo; figure w/ power vs. velo

% ----- DEREFERENCE BLOCK
v_vec            = linspace(V_stall,V_max,100);

weight_ttl       = UAVmin.weight.total;
weight_fuel      = UAVmin.fuel.W;
wing             = UAVmin.wing;
airfoilw         = UAVmin.airfoilw;
airfoilh         = UAVmin.airfoilh;
fuse             = UAVmin.fuse;
htail            = UAVmin.htail;
vtail            = UAVmin.vtail;
x_cg_emp         = UAVmin.stab.x_cg_empty;
x_cg_full        = UAVmin.stab.x_cg_full;
static_marg_emp  = UAVmin.stab.static_margin_empty;
static_marg_full = UAVmin.stab.static_margin_full;
rpm              = UAVmin.engn.rpm;
engn_HP          = UAVmin.engn.HP;
D                = UAVmin.prop.D;

DRAG.alt1000.empty = calc_drag_fn(v_vec,atmos(1).altitude,weight_ttl-weight_fuel,...
                wing,airfoilw,airfoilh,fuse,htail,vtail,x_cg_emp,static_marg_emp);
DRAG.alt1000.full = calc_drag_fn(v_vec,atmos(1).altitude,weight_ttl,...
                wing,airfoilw,airfoilh,fuse,htail,vtail,x_cg_full,static_marg_full);   
DRAG.alt7500.empty = calc_drag_fn(v_vec,atmos(2).altitude,weight_ttl-weight_fuel,...
                wing,airfoilw,airfoilh,fuse,htail,vtail,x_cg_emp,static_marg_emp);
DRAG.alt7500.full = calc_drag_fn(v_vec,atmos(2).altitude,weight_ttl,...
                wing,airfoilw,airfoilh,fuse,htail,vtail,x_cg_full,static_marg_full);

eta_p = calc_propeller(propeller,v_vec,rpm,D);
P_avail = engn_HP*eta_p;

% DRAG PROFILE
figure()
plot(v_vec, DRAG.alt1000.empty.D_t,'b'); hold on, grid on
plot(v_vec, DRAG.alt1000.full.D_t,'--b');
plot(v_vec, DRAG.alt7500.empty.D_t,'r');
plot(v_vec, DRAG.alt7500.full.D_t,'--r');
line([V_stall V_stall],[min(DRAG.alt1000.empty.D_t)-0.75 max(DRAG.alt7500.full.D_t)+0.75],'color','k');
line([V_loiter V_loiter ],[min(DRAG.alt1000.empty.D_t)-0.75 max(DRAG.alt7500.full.D_t)+0.75],'color','k');
line([V_cruise V_cruise],[min(DRAG.alt1000.empty.D_t)-0.75 max(DRAG.alt7500.full.D_t)+0.75],'color','k');
line([V_max   V_max  ],[min(DRAG.alt1000.empty.D_t)-0.75 max(DRAG.alt7500.full.D_t)+0.75],'color','k');
legend({'1000ft-empty','1000ft-full', '7500ft-empty', '7500ft-full'},'FontSize',10,'Location','SouthEast');
title('Drag Profile at Different Flight Conditions');   
xlabel('Velocity (ft/s)');  ylabel('Total Drag (lb_F)');
ylim([min(DRAG.alt1000.empty.D_t)-0.75 max(DRAG.alt7500.full.D_t)+0.75]);

% POWER PROFILE
figure()
plot(v_vec, DRAG.alt1000.empty.P_t,'b'); hold on, grid on
plot(v_vec, DRAG.alt1000.full.P_t,'--b');
plot(v_vec, DRAG.alt7500.empty.P_t,'r');
plot(v_vec, DRAG.alt7500.full.P_t,'--r');
plot(v_vec, P_avail,'k');
legend({'P_{reqd}:1000ft-empty','P_{reqd}:1000ft-full', 'P_{reqd}:7500ft-empty', 'P_{reqd}:7500ft-full','P_{avail}'},'FontSize',10,'Location','SouthEast');
title('Power Profile at Different Flight Conditions');
xlabel('Velocity (ft/s)');  ylabel('Total Power (HP)');

% LIFT TO DRAG
figure()
plot(DRAG.alt1000.empty.C_Dt,DRAG.alt1000.empty.C_L,'b'); hold on; grid on
plot(DRAG.alt1000.full.C_Dt,DRAG.alt1000.full.C_L,'--b');
plot(DRAG.alt7500.empty.C_Dt,DRAG.alt7500.empty.C_L,'r');
plot(DRAG.alt7500.full.C_Dt,DRAG.alt7500.full.C_L,'--r');
legend({'1000ft-empty','1000ft-full', '7500ft-empty', '7500ft-full'},'FontSize',10);
title('C_L vs. C_D Profile');
xlabel('C_D');  ylabel('C_L');
end

