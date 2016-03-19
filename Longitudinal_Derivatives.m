% Longitudinal Derivative

clear all; close all; clc;

S = 8.521;        % ft^2
b = 9.9967;       % ft
c = 0.8524;        % ft mean chord length
rho = 0.021;
g = 32.174;      % ft/s^2

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
m = W/g; % slug

% Stability Coefficient
C_Du = 0;       % no data assume 0
C_D0 = 0.0382 ; %0.060; %pioneer
C_Da = 0.0091; 

C_Lu = 0;      % no data assume 0
C_L0 = 0.0512; % 0.385; % pioneer

C_Ladot = 2.2892; %  2.42; % pioneer

C_mu = 0;      % no data assume 0
C_ma = -6.1293; % -2.12; % pioneer
C_mq =  -41.8653; %  -36.6; % pioneer
C_madot = -8.3731; % -11.0; % pioneer

% Dimensional Derivatives
X_u = -(C_Du+2*C_D0)*q*S/(m*u_0);
X_a = -(C_Da+C_L0)*q*S/(m);

Z_u = -(C_Lu+2*C_L0)*q*S/(m*u_0);
Z_a = -C_Ladot*c/(2*u_0)*q*S/m;

M_u = C_mu*(q*S*c)/(u_0*I_y);
M_a = C_ma*(q*S*c)/I_y;
M_q = C_mq*c/(2*u_0)*(q*S*c)/I_y;
M_adot = C_madot*c/(2*u_0)*(q*S*c)/I_y;


% Lateral-Directional Motion
lon_aero = [X_u                 X_a                 0          -g;
            Z_u/u_0             Z_a/u_0             1           0;
            M_u+M_adot*Z_u/u_0  M_a+M_adot*Z_a/u_0  M_q+M_adot  0;
            0                   0                   1           0];


lon_real_eig = eig(lon_aero);
real_roots_4th = real(lon_real_eig);
imag_roots_4th = imag(lon_real_eig);

% roll_approx = complex(L_p)
% sp_approx = complex((L_b*N_r - L_r*N_b)/L_b)
% 
% dutch_matrix = [Y_b/u_0   -(1);
%                 N_b       N_r ];
% dutch_approx = eig(dutch_matrix)

figure
hold on; grid on;
plot(lon_real_eig,'sr','MarkerFaceColor','y');
% plot(roll_approx,'sk','MarkerFaceColor','g');
% plot(sp_approx,'sk','MarkerFaceColor','g');
% plot(dutch_approx,'sk','MarkerFaceColor','g');
legend('Real Values','Location','NorthWest');
% whitebg('k');
grid on, 

damp(lon_real_eig)
% damp(roll_approx)
% damp(sp_approx)
% damp(dutch_approx)
% 
% freq_nat_dutch = sqrt((Y_b*N_r - N_b*0 + u_0*N_b)/u_0)
% damp_dutch = (-1/(2*freq_nat_dutch))*((Y_b+u_0*N_r)/u_0)