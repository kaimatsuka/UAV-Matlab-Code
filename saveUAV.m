function [UAV] = saveUAV(wing, airfoilw, fuse, htail, airfoilh, vtail, airfoilv,...
                    engn, fsys, prop, payld, stab, loadfact, WEIGHT, status)
% DESCRIPTION:
%   This funciton saves UAV parameters into UAV structure.
%
% INPUT:
%   wing = wing structure of current case
%   airfoilw = airfoil structure
%   fuse = fuselage structure
%   htail = horizontal tail structure
%   airfoilh = horizontal tail airfoil structure
%   vtail = vertical tail structure
%   airfoilv = vertical tail airfoil structure
%   engn = engine strucutre
%   fsys = fuel system structure
%   prop = propeller structure
%   payld = payload structure
%   stab = stability structure
%   loadfact    = load factor structure
%   WEIGHT = weight of aircraft
%
% REVISION HISTORY:
%   02/28: File created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    UAV.wing     = wing;
    UAV.airfoilw = airfoilw;
    UAV.fuse     = fuse;
    UAV.htail    = htail;
    UAV.airfoilh = airfoilh;
    UAV.vtail    = vtail;
    UAV.airfoilv = airfoilv;
    UAV.engn     = engn;
    UAV.fsys     = fsys;
    UAV.prop     = prop;
    UAV.payld    = payld;
    UAV.stab     = stab;
    UAV.loadfact = loadfact;
    UAV.weight   = WEIGHT; 
    UAV.status   = status;
end