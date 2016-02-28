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
            sfcl.x_cg_wing;   % ailron servo
            sfcl.x_cg_htail;  % rudder servo 
            sfcl.x_cg_vtail]; % elevator servo

xm_vec = x_cg_vec.*w_detail_vec; % (ft*lb)
stab.x_cg_full = sum(xm_vec)/sum(w_detail_vec); %(ft)


% X-direction - EMPTY FUEL
x_cg_empty_vec = [wing.x_cg; % wing CG
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
                sfcl.x_cg_wing;   % ailron servo
                sfcl.x_cg_htail;  % rudder servo 
                sfcl.x_cg_vtail]; % elevator servo
w_detail_emp_vec = w_detail_vec;
w_detail_emp_vec(8) = 0;    % Corresponds to fuel weight
xm_empty_vec = x_cg_empty_vec.*w_detail_emp_vec; %(ft*lb)
stab.x_cg_empty = sum(xm_empty_vec)/sum(w_detail_emp_vec); %(ft)

clearvars x_cg_empty_vec x_cg_vec xm_vec xm_empty_vec