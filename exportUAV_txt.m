function [ ] = exportUAV_txt( UAVstruct,filename,airfoils,engines)
%% exportUAV_txt.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:
%   This function exports aircraft parameters into a .txt file
%   
% INPUTS:
%   UAVstruct : structure of UAV to export

fileID = fopen(filename,'w');

%%%%% GEOMETRIC PROPERTIES OF UAV %%%%%
%--- WING
fprintf(fileID,'-----------GEOMETRIC PROPERTIES----------\n');
fprintf(fileID,'-----WING PROPERTIES-----\n');
fprintf(fileID,'Airfoil = %8s\n',char(airfoils(UAVstruct.airfoilw.ind).name));
fprintf(fileID,'Wing Area (S) = %4.4f ft\n',UAVstruct.wing.S);
fprintf(fileID,'Aspect Ratio (A) = %4.4f ft\n',UAVstruct.wing.A);
fprintf(fileID,'Taper Ratio (lam) = %4.4f\n',UAVstruct.wing.lam);
fprintf(fileID,'Wing Quarter Chord Sweep (deg) (lam_q) = %4.4f\n',UAVstruct.wing.lam_q);
fprintf(fileID,'Distance from Nose to Quarter Chord Pt (h_q) = %4.4f\n',UAVstruct.wing.h_q);
fprintf(fileID,'Distance from Nose to Root Leading Edge (h_LE) = %4.4f ft\n',UAVstruct.wing.x_LE);
fprintf(fileID,'Oswald Efficiency (e) = %4.4f\n',UAVstruct.wing.e);
fprintf(fileID,'Wing Span (b) = %4.4f ft\n',UAVstruct.wing.b);
fprintf(fileID,'Chord Length (c) = %4.4f ft\n',UAVstruct.wing.c);
fprintf(fileID,'Thickness (t) = %4.4f ft\n',UAVstruct.wing.t);
fprintf(fileID,'Maximum Thickness Ratio (t/c) = %4.4f\n', UAVstruct.wing.mtr);
fprintf(fileID,'Maximum Thickness Location (mtl) = %4.4f\n',UAVstruct.wing.mtl);
fprintf(fileID,'Wing Root Chord Length (c_r) = %4.4f ft\n',UAVstruct.wing.c_r);
fprintf(fileID,'Wing Tip Chord Length (c_t) = %4.4f ft\n',UAVstruct.wing.c_t);
fprintf(fileID,'Max Thickness at Root (t_r) = %4.4f ft\n',UAVstruct.wing.t_r);
fprintf(fileID,'Max Thickness at Tip (t_t) = %4.4f ft\n',UAVstruct.wing.t_t);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft\n',UAVstruct.wing.x_cg);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.wing.z_cg);
fprintf(fileID,'CL_max = %4.4f\n',UAVstruct.airfoilw.CLmax);
fprintf(fileID,'Lift Curve Slope (a_w) = %4.4f (1/rad)\n',UAVstruct.airfoilw.a_w);

%--- FUSELAGE
fprintf(fileID,'\n-----FUSELAGE PROPERTIES-----\n');
fprintf(fileID,'Fuselage Length (L) = %4.4f ft\n',UAVstruct.fuse.L);
fprintf(fileID,'Fuselage Diameter (D) = %4.4f ft\n',UAVstruct.fuse.D);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft\n',UAVstruct.fuse.x_cg);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft\n',UAVstruct.fuse.z_cg);

%--- HORIZONTAL TAIL
fprintf(fileID,'\n-----HORIZONTAL TAIL PROPERTIES-----\n');
fprintf(fileID,'Airfoil = %8s\n',char(airfoils(UAVstruct.airfoilh.ind).name));
fprintf(fileID,'Horizontal Tail Area (S) = %4.4f ft\n',UAVstruct.htail.S);
fprintf(fileID,'Aspect Ratio (A) = %4.4f ft\n',UAVstruct.htail.A);
fprintf(fileID,'Taper Ratio (lam) = %4.4f\n',UAVstruct.htail.lam);
fprintf(fileID,'Horizontal Tail Quarter Chord Sweep (deg) (lam_q) = %4.4f\n',UAVstruct.htail.lam_q);
fprintf(fileID,'Distance from Nose to Quarter Chord Pt (h) = %4.4f\n',UAVstruct.htail.h);
fprintf(fileID,'Distance from Nose to Root Leading Edge (h_LE) = %4.4f ft\n',UAVstruct.htail.x_LE);
fprintf(fileID,'Oswald Efficiency (e) = %4.4f\n',UAVstruct.htail.e);
fprintf(fileID,'Horizontal Tail Span (b) = %4.4f ft\n',UAVstruct.htail.b);
fprintf(fileID,'Chord Length (c) = %4.4f ft\n',UAVstruct.htail.c);
fprintf(fileID,'Thickness (t) = %4.4f ft\n',UAVstruct.htail.mtr*UAVstruct.htail.c);
fprintf(fileID,'Maximum Thickness Ratio (t/c) = %4.4f\n', UAVstruct.htail.mtr);
fprintf(fileID,'Maximum Thickness Location (mtl) = %4.4f\n',UAVstruct.htail.mtl);
fprintf(fileID,'Horizontal Tail Root Chord Length (c_r) = %4.4f ft\n',UAVstruct.htail.c_r);
fprintf(fileID,'Horizontal Tail Tip Chord Length (c_t) = %4.4f ft\n',UAVstruct.htail.c_t);
fprintf(fileID,'Max Thickness at Root (t_r) = %4.4f ft\n',UAVstruct.htail.t_r);
fprintf(fileID,'Max Thickness at Tip (t_t) = %4.4f ft\n',UAVstruct.htail.t_t);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft\n',UAVstruct.htail.x_cg);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.wing.z_cg);
fprintf(fileID,'CL_max = %4.4f\n',UAVstruct.airfoilw.CLmax);
fprintf(fileID,'Lift Curve Slope (a_t) = %4.4f (1/rad)\n',UAVstruct.airfoilh.a_t);

%--- VERTICAL TAIL
fprintf(fileID,'\n-----VERTICAL TAIL PROPERTIES-----\n');
fprintf(fileID,'Airfoil = %8s\n',char(airfoils(UAVstruct.airfoilv.ind).name));
fprintf(fileID,'Horizontal Tail Area (S) = %4.4f ft\n',UAVstruct.vtail.S);
fprintf(fileID,'Aspect Ratio (A) = %4.4f ft\n',UAVstruct.vtail.A);
fprintf(fileID,'Taper Ratio (lam) = %4.4f\n',UAVstruct.vtail.lam);
fprintf(fileID,'Horizontal Tail Quarter Chord Sweep (deg) (lam_q) = %4.4f\n',UAVstruct.vtail.lam_q);
fprintf(fileID,'Distance from Nose to Quarter Chord Pt (h) = %4.4f\n',UAVstruct.vtail.h);
fprintf(fileID,'Distance from Nose to Root Leading Edge (h_LE) = %4.4f ft\n',UAVstruct.vtail.x_LE);
fprintf(fileID,'Oswald Efficiency (e) = %4.4f\n',UAVstruct.vtail.e);
fprintf(fileID,'Horizontal Tail Span (b) = %4.4f ft\n',UAVstruct.vtail.b);
fprintf(fileID,'Chord Length (c) = %4.4f ft\n',UAVstruct.vtail.c);
fprintf(fileID,'Thickness (t) = %4.4f ft\n',UAVstruct.vtail.t);
fprintf(fileID,'Maximum Thickness Ratio (t/c) = %4.4f\n', UAVstruct.vtail.mtr);
fprintf(fileID,'Maximum Thickness Location (mtl) = %4.4f\n',UAVstruct.vtail.mtl);
fprintf(fileID,'Horizontal Tail Root Chord Length (c_r) = %4.4f ft\n',UAVstruct.vtail.c_r);
fprintf(fileID,'Horizontal Tail Tip Chord Length (c_t) = %4.4f ft\n',UAVstruct.vtail.c_t);
fprintf(fileID,'Max Thickness at Root (t_r) = %4.4f ft\n',UAVstruct.vtail.t_r);
fprintf(fileID,'Max Thickness at Tip (t_t) = %4.4f ft\n',UAVstruct.vtail.t_t);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft\n',UAVstruct.vtail.x_cg);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.vtail.z_cg);
fprintf(fileID,'CL_max = %4.4f\n',UAVstruct.airfoilv.CLmax);
fprintf(fileID,'Lift Curve Slope (a_t) = %4.4f (1/rad)\n',UAVstruct.airfoilv.a_v);

%---Engine
fprintf(fileID,'\n-----ENGINE PROPERTIES-----\n');
fprintf(fileID,'Engine Name = %10s\n',char(engines(UAVstruct.engn.ind).name));
fprintf(fileID,'Power (HP) = %4.4f\n',UAVstruct.engn.HP);
fprintf(fileID,'RPM = %4.4f \n',UAVstruct.engn.rpm);
fprintf(fileID,'Engine Weight (lb) = %4.4f \n',UAVstruct.engn.W_bare);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft \n',UAVstruct.engn.x_cg);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.engn.z_cg);

%---Propeller
fprintf(fileID,'\n-----PROPELLER PROPERTIES-----\n');
fprintf(fileID,'Diameter = %4.4f ft',UAVstruct.prop.D);
fprintf(fileID,'Pitch = %4.4f deg',UAVstruct.prop.pitch);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft \n',UAVstruct.prop.x_cg);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.prop.z_cg);

%---Fuel % Fuel System
fprintf(fileID,'\n-----FUEL & FUEL SYSTEM-----\n');
fprintf(fileID,'Fuel Weight = %4.4f lb\n',UAVstruct.fuel.W);
fprintf(fileID,'Fuel System Weight = %4.4f lb\n',UAVstruct.fsys.W);
fprintf(fileID,'Fuel System Diameter = %4.4f ft\n',UAVstruct.fsys.D);
fprintf(fileID,'Fuel System Length = %4.4f ft\n',UAVstruct.fsys.length);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft \n',UAVstruct.fsys.x_cg);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.fsys.z_cg);

%---Payload Specs
fprintf(fileID,'\n-----PAYLOAD SPECS-----\n');
fprintf(fileID,'Total Weight = %4.4f lbs\n',UAVstruct.payld.w_total);
fprintf(fileID,'--EO/IR Sensor\n');
fprintf(fileID,'Weight = %4.4f lbs\n',UAVstruct.payld.w_EOIR);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft \n',UAVstruct.payld.x_cg_EOIR);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.payld.z_cg_EOIR);
fprintf(fileID,'--LiDAR Sensor\n');
fprintf(fileID,'Weight = %4.4f lbs\n',UAVstruct.payld.w_LiDAR);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft \n',UAVstruct.payld.x_cg_LiDAR);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.payld.z_cg_LiDAR);
fprintf(fileID,'--SAR Sensor\n');
fprintf(fileID,'Weight = %4.4f lbs\n',UAVstruct.payld.w_SAR);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft \n',UAVstruct.payld.x_cg_SAR);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.payld.z_cg_SAR);
fprintf(fileID,'--Wave Relay\n');
fprintf(fileID,'Weight = %4.4f lbs\n',UAVstruct.payld.w_WR);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft \n',UAVstruct.payld.x_cg_WR);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.payld.z_cg_WR);
fprintf(fileID,'--IMU\n');
fprintf(fileID,'Weight = %4.4f lbs\n',UAVstruct.payld.w_IMU);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft \n',UAVstruct.payld.x_cg_IMU);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.payld.z_cg_IMU);
fprintf(fileID,'--SAR Sensor\n');
fprintf(fileID,'Weight = %4.4f lbs\n',UAVstruct.payld.w_SAR);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft \n',UAVstruct.payld.x_cg_SAR);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.payld.z_cg_SAR);
fprintf(fileID,'--Antenna\n');
fprintf(fileID,'Weight = %4.4f lbs\n',UAVstruct.payld.w_ANT);
fprintf(fileID,'X-Center of Gravity (x_cg) = %4.4f ft \n',UAVstruct.payld.x_cg_ANT);
fprintf(fileID,'Z-Center of Gravity (z_cg) = %4.4f ft \n',UAVstruct.payld.z_cg_ANT);

fprintf(fileID,'\n-----------WEIGHT PROPERTIES----------\n');
fprintf(fileID,'Total Weight = %4.4f lbs\n',UAVstruct.weight.total);
fprintf(fileID,'  Wing = %4.4f lbs\n',UAVstruct.weight.wing);
fprintf(fileID,'  Fuselage = %4.4f lbs\n',UAVstruct.weight.fuse);
fprintf(fileID,'  Horizontal Tail = %4.4f lbs\n',UAVstruct.weight.htail);
fprintf(fileID,'  Vertical Tail = %4.4f lbs\n',UAVstruct.weight.vtail);
fprintf(fileID,'  Fuel System = %4.4f lbs\n',UAVstruct.weight.fsys);
fprintf(fileID,'  Engine = %4.4f lbs\n',UAVstruct.weight.engn);
fprintf(fileID,'  Avionics/Payload = %4.4f lbs\n',UAVstruct.weight.avion);
fprintf(fileID,'  Fuel = %4.4f lbs\n',UAVstruct.weight.fuel);
fprintf(fileID,'  Propeller = %4.4f lbs\n',UAVstruct.weight.prop);

fprintf(fileID,'\n-----------PERFORMANCE----------\n');
fprintf(fileID,'\n-----------LIFT/DRAG/POWER PROPERTIES----------\n');
fprintf(fileID,'--Condition: Altitude 1000ft, Fuel: Full\n');
fprintf(fileID,'Speed (ft/s):       %10.4f (Stall)   %10.4f (Loiter)\n',UAVstruct.drag.alt1000.full.v);
fprintf(fileID,'Lift Coeff (C_L):   %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.C_L);
fprintf(fileID,'  C_Lw:             %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.C_Lw);
fprintf(fileID,'  C_Lh:             %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.C_Lh);
fprintf(fileID,'Drag Coeff (C_D):   %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.C_Dt);
fprintf(fileID,'  C_Dp:             %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.C_Dp);
fprintf(fileID,'  C_Diw:            %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.C_Diw);
fprintf(fileID,'  C_Dih:            %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.C_Dih);
fprintf(fileID,'  C_Dairf:          %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.C_Dairf);
fprintf(fileID,'  C_Dtrim:          %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.C_Dt-UAVstruct.drag.alt1000.full.C_Dp-UAVstruct.drag.alt1000.full.C_Diw-UAVstruct.drag.alt1000.full.C_Dairf);
fprintf(fileID,'Moment AC (Cm_ac):  %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.Cm_ac);
fprintf(fileID,'Trim Angle (i_t):   %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.i_t);
fprintf(fileID,'Total Drag (lbf):   %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.D_t);
fprintf(fileID,'  Parasite (lbf):   %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.D_p);
fprintf(fileID,'  Induced (lbf):    %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.D_iw);
fprintf(fileID,'  Airfoil (lbf):    %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.D_airf);
fprintf(fileID,'  Trim (lbf):       %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.D_t-UAVstruct.drag.alt1000.full.D_p-UAVstruct.drag.alt1000.full.D_iw-UAVstruct.drag.alt1000.full.D_airf);
fprintf(fileID,'Power Reqd (HP):    %10.4f           %10.4f\n',UAVstruct.drag.alt1000.full.P_t);

fprintf(fileID,'\n--Condition: Altitude 1000ft, Fuel: Empty\n');
fprintf(fileID,'Speed (ft/s):       %10.4f (Stall)   %10.4f (Loiter)\n',UAVstruct.drag.alt1000.empty.v);
fprintf(fileID,'Lift Coeff (C_L):   %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.C_L);
fprintf(fileID,'  C_Lw:             %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.C_Lw);
fprintf(fileID,'  C_Lh:             %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.C_Lh);
fprintf(fileID,'Drag Coeff (C_D):   %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.C_Dt);
fprintf(fileID,'  C_Dp:             %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.C_Dp);
fprintf(fileID,'  C_Di:             %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.C_Di);
fprintf(fileID,'  C_Dairf:          %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.C_Dairf);
fprintf(fileID,'  C_Dtrim:          %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.C_Dt-UAVstruct.drag.alt1000.empty.C_Dp-UAVstruct.drag.alt1000.empty.C_Di-UAVstruct.drag.alt1000.empty.C_Dairf);
fprintf(fileID,'Moment AC (Cm_ac):  %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.Cm_ac);
fprintf(fileID,'Trim Angle (i_t):   %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.i_t);
fprintf(fileID,'Total Drag (lbf):   %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.D_t);
fprintf(fileID,'  Parasite (lbf):   %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.D_p);
fprintf(fileID,'  Induced (lbf):    %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.D_i);
fprintf(fileID,'  Airfoil (lbf):    %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.D_airf);
fprintf(fileID,'  Trim (lbf):       %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.D_t-UAVstruct.drag.alt1000.empty.D_p-UAVstruct.drag.alt1000.empty.D_i-UAVstruct.drag.alt1000.empty.D_airf);
fprintf(fileID,'Power Reqd (HP):    %10.4f           %10.4f\n',UAVstruct.drag.alt1000.empty.P_t);

fprintf(fileID,'\n--Condition: Altitude 7500ft, Fuel: Full\n');
fprintf(fileID,'Speed (ft/s):       %10.4f (Stall)   %10.4f (Loiter)\n',UAVstruct.drag.alt7500.full.v);
fprintf(fileID,'Lift Coeff (C_L):   %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.C_L);
fprintf(fileID,'  C_Lw:             %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.C_Lw);
fprintf(fileID,'  C_Lh:             %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.C_Lh);
fprintf(fileID,'Drag Coeff (C_D):   %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.C_Dt);
fprintf(fileID,'  C_Dp:             %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.C_Dp);
fprintf(fileID,'  C_Di:             %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.C_Di);
fprintf(fileID,'  C_Dairf:          %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.C_Dairf);
fprintf(fileID,'  C_Dtrim:          %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.C_Dt-UAVstruct.drag.alt7500.full.C_Dp-UAVstruct.drag.alt7500.full.C_Di-UAVstruct.drag.alt7500.full.C_Dairf);
fprintf(fileID,'Moment AC (Cm_ac):  %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.Cm_ac);
fprintf(fileID,'Trim Angle (i_t):   %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.i_t);
fprintf(fileID,'Total Drag (lbf):   %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.D_t);
fprintf(fileID,'  Parasite (lbf):   %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.D_p);
fprintf(fileID,'  Induced (lbf):    %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.D_i);
fprintf(fileID,'  Airfoil (lbf):    %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.D_airf);
fprintf(fileID,'  Trim (lbf):       %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.D_t-UAVstruct.drag.alt7500.full.D_p-UAVstruct.drag.alt7500.full.D_i-UAVstruct.drag.alt7500.full.D_airf);
fprintf(fileID,'Power Reqd (HP):    %10.4f           %10.4f\n',UAVstruct.drag.alt7500.full.P_t);

fprintf(fileID,'\n--Condition: Altitude 7500ft, Fuel: Empty\n');
fprintf(fileID,'Speed (ft/s):       %10.4f (Stall)   %10.4f (Loiter)\n',UAVstruct.drag.alt7500.empty.v);
fprintf(fileID,'Lift Coeff (C_L):   %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.C_L);
fprintf(fileID,'  C_Lw:             %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.C_Lw);
fprintf(fileID,'  C_Lh:             %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.C_Lh);
fprintf(fileID,'Drag Coeff (C_D):   %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.C_Dt);
fprintf(fileID,'  C_Dp:             %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.C_Dp);
fprintf(fileID,'  C_Di:             %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.C_Di);
fprintf(fileID,'  C_Dairf:          %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.C_Dairf);
fprintf(fileID,'  C_Dtrim:          %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.C_Dt-UAVstruct.drag.alt7500.empty.C_Dp-UAVstruct.drag.alt7500.empty.C_Di-UAVstruct.drag.alt7500.empty.C_Dairf);
fprintf(fileID,'Moment AC (Cm_ac):  %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.Cm_ac);
fprintf(fileID,'Trim Angle (i_t):   %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.i_t);
fprintf(fileID,'Total Drag (lbf):   %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.D_t);
fprintf(fileID,'  Parasite (lbf):   %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.D_p);
fprintf(fileID,'  Induced (lbf):    %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.D_i);
fprintf(fileID,'  Airfoil (lbf):    %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.D_airf);
fprintf(fileID,'  Trim (lbf):       %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.D_t-UAVstruct.drag.alt7500.empty.D_p-UAVstruct.drag.alt7500.empty.D_i-UAVstruct.drag.alt7500.empty.D_airf);
fprintf(fileID,'Power Reqd (HP):    %10.4f           %10.4f\n',UAVstruct.drag.alt7500.empty.P_t);

fprintf(fileID,'\n-----------PROPELLER EFFICIENCY----------\n');

end

