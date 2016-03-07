function [eta_p] = calc_propeller(propeller,V,RPM,D)
% calc_propeller.m
%
% DESCRITPION:
%   calculate propeller efficiency
%
% INPUTS:
%   propeller.J     = array of v/nD 1xN
%            .eta_p = array of propeller efficiency corresponding to J 1xN
%   V_air    = UAV velocity (ft/s) 1xM
%   RPM      = engine rpm (rpm) 1x1
%   engn.HP  = engine horse power (HP) 1x1
%
% OUTPUTS:
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_this = V/(RPM/60*D);
eta_p = interp1(propeller.J,propeller.eta_p,J_this);


