function out = get_coeffs(in)

alpha = in(1);
alpha_dot = in(2);
beta = in(3);
elevator = in(4);
aileron = in(5);
p = in(6);
q = in(7);
r = in(8);
rudder = in(9);
Vt = in(10);

cbar =  1.80;     % [ft]        Bray pg 31
bw   =  16.90;    % [ft]        Bray pg 31


% aero data - linear stability coefficients
CL0     =    SDERIV.CL0;         % []          Bray pg 33
CL_a    =    SDERIV.CL_a;        % [/rad]      Bray pg 33
CL_adot =    SDREIV.CL_adot;     % [/rad]      Bray pg 33
CL_q    =    SDERIV.CL_q;        % [/rad]      Bray pg 33
CL_de   =    SDERIV.CL_de;       % [/rad]      Bray pg 33
%CLfa CLfa.dat 0 1               %             Bray pg 50, Table 4.7
%CLfade CLfade.dat 0 1 1         %             Bray pg 41, Table 4.4

CD0     =    SDERIV.CD0;         % []          Bray pg 33
CD_a    =    SDERIV.CD_a;        % [/rad]      Bray pg 33
CD_de   =    SDERIV.CD_de;       % [/rad]      Bray pg 33
%CDfa Cdfa.dat 0 1               %             Bray pg 50, Table 4.7
%CDfade CDfade.dat 0 1 1         %             Bray pg 45, Table 4.6

Cm0     =   SDERIV.Cm0;          % []          Bray pg 33
Cm_a    =   SDERIV.Cm_a;         % [/rad]      Bray pg 33
Cm_adot =   SDERIV.Cm_adot;      % [/rad]      Bray pg 33
Cm_q    =   SDERIV.Cm_q;         % [/rad]      Bray pg 33
Cm_de   =   SDERIV.Cm_de;        % [/rad]      Bray pg 33
%Cmfade Cmfade.dat 0 1 1         %             Bray pg 34, Figure 4.2

CY_beta =   SDERIV.CY_beta;      % [/rad]      Bray pg 33
%CY_p      <CY_p>                % [/rad]      no data
%CY_r      <CY_r>                % [/rad]      no data
%CY_da     <CY_da>               % [/rad]      no data
CY_dr   =   SDERIV.CY_dr;        % [/rad]      Bray pg 33
%CYfada <CYfada.dat>             % []          no data
%CYfbetadr CYfbetadr.dat 0 1 1   % []          Bray pg 62, Figure 4.19

Cl_beta =   SDERIV.Cl_beta;      % [/rad]      Bray pg 33
Cl_p    =   SDERIV.Cl_p;         % [/rad]      Bray pg 33
Cl_r    =   SDERIV.Cl_r;         % [/rad]      Bray pg 33
Cl_da   =   SDERIV.Cl_da;        % [/rad]      Bray pg 33 (sign reversed)
Cl_dr   =   SDERIV.Cl_dr;        % [/rad]      Bray pg 33
%Clfada Clfada.dat 0 1 1         %             Bray pg 58, Table 4.8

Cn_beta =   SDERIV.Cn_beta;      % [/rad]      Bray pg 33
Cn_p    =   SDERIV.Cn_p;         % [/rad]      Bray pg 33
Cn_r    =   SDERIV.Cn_r;         % [/rad]      Bray pg 33
Cn_da   =   SDERIV.Cn_da;        % [/rad]      Bray pg 33 (sign reversed)
Cn_dr   =   SDERIV.Cn_dr;        % [/rad]      Bray pg 33
%Cnfada Cnfada.dat 0 1 1         %             Bray pg 61, Table 4.9
%Cnfbetadr Cnfbetadr.dat 0 1 1   %             Bray pg 63, Figure 4.20


% Calculate CL, CD, CM

CL = CL0 + CL_a*alpha + CL_de*elevator +...
     CL_adot*alpha_dot*cbar/(2*Vt)+CL_q*q*cbar/(2*Vt);
 
CM = Cm0 + Cm_a*alpha + Cm_de*elevator +... 
     Cm_adot*alpha_dot*cbar/(2*Vt)+Cm_q*q*cbar/(2*Vt);
 
CD = CD0 + CD_a*alpha + CD_de*elevator;

% Calculate CY, Cl, CN

CY = CY_beta*beta + CY_dr*rudder;

Cl = Cl_beta*beta + Cl_da*aileron + Cl_dr*rudder + ...
    Cl_p*p*bw/(2*Vt)+Cl_r*r*bw/(2*Vt);

CN = Cn_beta*beta + Cn_da*aileron + Cn_dr*rudder + ...
    Cn_p*p*bw/(2*Vt)+Cn_r*r*bw/(2*Vt);

out(1) = CD;
out(2) = CY;
out(3) = CL;
out(4) = Cl;
out(5) = CM;
out(6) = CN;

 