function [static_margin] = calc_static_marg (x_cg, x_n, wing)
% DESCRIPTION:
%   This file determines static margin at full fuel and empty fuel
%
% INPUT:
%   x_cg: center of gravity relative to nose (ft)
%   x_n: neutral point relative to nose (ft)
%   wing: wing structure
%
% OUTPUT:
%   static_margin: static margin @ full fuel
%
% REVISION HISTORY:
%   02/27: File created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_cg_LE  = x_cg - (wing.x_cg - (0.25*wing.c)); % Distance from wing LE to CG (ft)
    h_cg_LE  = x_cg_LE/wing.c;  % Normalized distance from LE to CG [rel to chord] (-) 
    
    x_np_LE   = x_n - (wing.x_cg - (0.25*wing.c));  % Distance from wing LE to NP (ft)
    h_np_LE   = x_np_LE/wing.c; % Normalized distance from LE to NP [rel to chord] (-)
    
    static_margin = h_np_LE - h_cg_LE;
end