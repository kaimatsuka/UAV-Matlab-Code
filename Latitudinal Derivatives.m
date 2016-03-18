clear all; close all; clc;
S = 8.521;        % ft^2
b = 9.9967;       % ft
c = 0.8524;        % ft mean chors
rho = 0.021;
g = 32.174;

% Normal Flight Condition
h = 7500;          % ft 
v_to = 75.9;     % ft/s
W = 32.8393;       % lb
h_cg = 0.295*c;   % 29.5% mean aero chord
I_x = 2.91;      % slug-ft^2 found from xflr5
I_y = 10.98;     % slug-ft^2 found from xflr5
I_z = 13.55;     % slug-ft^2 found from xflr5
I_xz = -0.455;   % slug-ft^2 found from xflr5

u_0 = 54.2669;
q = 0.5*rho*u_0^2;
m = W/g;

% Stability Coefficient
C_Yb    = -0.5939;    % y-force change due to side slip angle
C_Lp    = -0.9051;    % roll moment change due to roll rate
C_Lr    = 0.1254;    % roll moment change due to yaw rate
C_Lda   = 0.7636;    % roll moment change due to aileron deflection
C_Ydr   = 0.4751;    % y-force change due to rudder deflection
C_Lb    = 0;    % roll moment change due to side slip angle
C_Np    = -0.11;   % yaw moment change due to roll rate
C_Nr    = -0.0682;    % yaw moment change due to yaw rate
C_Nda   = 0.0036;   % yaw moment change due to aileron deflection
C_Ldr   = 0.0511;    % roll moment change due to rudder deflection
C_Nb    = 0.161;    % yaw moment change due to side slip angle
C_Ndr   = -0.161;    % yaw moment change due to rudder deflection

% Dimensional Derivatives
Y_b  = q*(S/m)*C_Yb;                  % ft/s^2
Y_p  = 0;                             % ft/s
Y_r  = 0;
Y_da = 0;
Y_dr = ((q*S)/m)*C_Ydr;               % ft/s^2

L_b  = ((q*S*b)/I_x)*C_Lb;            % 1/s^2
L_p  = ((q*S*b^2)/(2*I_x*u_0))*C_Lp;  % 1/s
L_r  = ((q*S*b^2)/(2*I_x*u_0))*C_Lr;  % 1/s
L_da = ((q*S*b)/I_x)*C_Lda;           % 1/s^2
L_dr = ((q*S*b)/I_x)*C_Ldr;           % 1/s^2

N_b  = ((q*S*b)/I_z)*C_Nb;            % 1/s^2
N_p  = ((q*S*b^2)/(2*I_z*u_0))*C_Np;  % 1/s
N_r  = ((q*S*b^2)/(2*I_z*u_0))*C_Nr;  % 1/s
N_da = ((q*S*b)/I_z)*C_Nda;           % 1/s^2
N_dr = ((q*S*b)/I_z)*C_Ndr;           % 1/s^2

% Lateral-Directional Motion
lat_aero = [Y_b/u_0  0      -(1)           g/u_0;
           L_b       L_p    L_r            0    ;
           N_b       N_p    N_r            0    ;
           0         1      0              0    ];
lat_cont = [0           Y_dr/u_0;
            L_da        L_dr    ;
            N_da        N_dr    ;
            0           0       ];

lat_real_eig = eig(lat_aero)
real_roots_4th = real(lat_real_eig);
imag_roots_4th = imag(lat_real_eig);

roll_approx = complex(L_p)
sp_approx = complex((L_b*N_r - L_r*N_b)/L_b)

dutch_matrix = [Y_b/u_0   -(1);
                N_b       N_r ];
dutch_approx = eig(dutch_matrix)

figure
hold on; grid on;
plot(lat_real_eig,'sr','MarkerFaceColor','y');
plot(roll_approx,'sk','MarkerFaceColor','g');
plot(sp_approx,'sk','MarkerFaceColor','g');
plot(dutch_approx,'sk','MarkerFaceColor','g');
legend('Real Values','Approximations','Location','NorthWest');
whitebg('black');

damp(lat_real_eig)
damp(roll_approx)
damp(sp_approx)
damp(dutch_approx)

freq_nat_dutch = sqrt((Y_b*N_r - N_b*0 + u_0*N_b)/u_0)
damp_dutch = (-1/(2*freq_nat_dutch))*((Y_b+u_0*N_r)/u_0)