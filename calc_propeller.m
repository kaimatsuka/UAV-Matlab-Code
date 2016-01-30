% calc_propeller.m
%
% DESCRITPION:
%   calculate propeller efficiency
%
% INPUTS:
%   engn.HP = engine horse power (HP)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prop.n_blades = 2;     %number of propeller blades
prop.N_prop = 7500/60; % rotational speed (selected from engine) (RPS)

if prop.n_blades==2
    prop.d_prop = 1.7*(engn.HP)^.25; % Raymer (pg 315)
elseif n_blades == 3
    prop.d_prop = 1.6*(engn.HP)^.25; % Raymer (pg 315)
elseif n_blades == 4
    prop.d_prop = 1.5*(engn.HP)^.25; % Raymer (pg 315)
end

% propeller efficiency 
%{
    Ap     = pi*d_prop^2/4; %[ft^2] disc area 
    T      = 3.5;  %[lbf] thrust force (= drag at this velocity)
    v_prop = [50:10:1000]; % specify velocity range for plotting

    k = 20; % ratio between Fy and Fx     http://www.rcex.cz/?p=3593
    x = v_prop/(pi*d_prop*N_prop);
    eta_p = (k-x)./(k*x+1).*x;
%}
%{
dV = sqrt(v_prop.^2+2*T/rho/Ap)-v_prop; %[ft/s] Gundlach (pg 292)
eta_p_ideal = 1./(dV/2./v_prop+1); 
eta_p_push = (1-0.12);           % push prop efficiency factor (Toohey)
eta_p_nonideal = 0.8; % non ideal efficiency 
eta_p = eta_p_ideal*eta_p_nonideal*eta_p_push;
%}

prop.eta_p = 0.8*(1-0.12); % propulsive efficiency (-12% for push prop: Toohey)