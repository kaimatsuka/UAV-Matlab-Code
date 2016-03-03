% calc_mass_properties.m
%
% DESCRIPTION:
%   
%
% INPUT:
%   
%
% OUTPUT:
%   stab.x_cg_full: center of gravity relative to nose - FULL FUEL (ft)
%   stab.x_cg_empty : CG relative to nose - EMPTY FUEL (ft)
%
% REVISION HISTORY:
%   02/16: File created.
%   02/27: Added stab structure; added empty fuel CG location
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coordinate System
%   Origin of coordinate system is located at the nose
%   X-direction: from nose to tail
%   Y-direction: from left to right
%   Z-direction: from bottom to top

% X-direction - FULL FUEL

x_cg_vec = [wing.x_cg; % wing CG
            fuse.x_cg;
            htail.x_cg;
            vtail.x_cg;
            engn.x_cg;
            fsys.x_cg;
            prop.x_cg;
            fsys.x_cg; % fuel CG ASSUME SAME CG AS FUEL SYSTEM
            payld.x_cg_EOIR;  % electro optical
            payld.x_cg_SAR;   % synthetic apateur 
            payld.x_cg_LiDAR; % LiDAR
            payld.x_cg_ANT;
            payld.x_cg_IMU;
            payld.x_cg_WR;
            sfcl.ail.x_cg;   % ailron servo
            sfcl.rudd.x_cg;  % rudder servo 
            sfcl.elev.x_cg]; % elevator servo

xm_vec = x_cg_vec.*w_detail_vec; % (ft*lb)
stab.x_cg_full = sum(xm_vec)/sum(w_detail_vec); %(ft)


% X-direction - EMPTY FUEL
w_detail_emp_vec = w_detail_vec;
w_detail_emp_vec(8) = 0;    % Corresponds to fuel weight
xm_empty_vec = x_cg_vec.*w_detail_emp_vec; %(ft*lb)
stab.x_cg_empty = sum(xm_empty_vec)/sum(w_detail_emp_vec); %(ft)

% Z-direction - FULL FUEL
z_cg_vec = [wing.z_cg; % wing CG
            fuse.z_cg;
            htail.z_cg;
            vtail.z_cg;
            engn.z_cg;
            fsys.z_cg;
            prop.z_cg;
            fsys.z_cg; % fuel CG ASSUME SAME CG AS FUEL SYSTEM
            payld.z_cg_EOIR;  % electro optical
            payld.z_cg_SAR;   % synthetic apateur 
            payld.z_cg_LiDAR; % LiDAR
            payld.z_cg_ANT;
            payld.z_cg_IMU;
            payld.z_cg_WR;
            sfcl.ail.z_cg;   % ailron servo
            sfcl.rudd.z_cg;  % rudder servo 
            sfcl.elev.z_cg]; % elevator servo

zm_vec = z_cg_vec.*w_detail_vec; % (ft*lb)
stab.z_cg_full = sum(zm_vec)/sum(w_detail_vec); %(ft)

w_detail_emp_vec = w_detail_vec;
w_detail_emp_vec(8) = 0;    % Corresponds to fuel weight
zm_empty_vec = z_cg_vec.*w_detail_emp_vec; %(ft*lb)
stab.z_cg_empty = sum(zm_empty_vec)/sum(w_detail_emp_vec); %(ft)

clearvars x_cg_empty_vec x_cg_vec xm_vec xm_empty_vec zm_vec zm_empty_vec