% calc_mass_properties.m
%
% DESCRIPTION:
%   
%
% INPUT:
%   
%
% OUTPUT:
%   
%
% REVISION HISTORY:
%   02/16: File created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coordinate System
%   Origin of coordinate system is located at the nose
%   X-direction: from nose to tail
%   Y-direction: from left to right
%   Z-direction: from bottom to top

% X-direction

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
            ail.x_cg_sfcl;   % ailron servo
            rudd.x_cg_sfcl;  % rudder servo 
            elev.x_cg_sfcl]; % elevator servo

xm_vec = x_cg_vec.*w_detail_vec; % (ft*lb)
x_cg_total = sum(xm_vec)/sum(w_detail_vec); %(ft)



