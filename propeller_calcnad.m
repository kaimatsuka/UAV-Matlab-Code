%*************************************************************************
% Numerical Interpolation
%*************************************************************************

clear all; close all; clc;

% Known parameters
n = 7500/60; % revPS
D = 2; % Diameter (ft)

% Reading in excel files
[num,txt] = xlsread('Propeller_Design.xlsx','ClarkY Efficiency', 'A1:B35','basic');

% Separating excel data, performing pchip interpolation
x = num(:, 1);
y15 = num(:, 2);
xx = num(1,1):0.001:num(length(num), 1);
yy = pchip(x,y15,xx);

% Plotting parameter vs efficiency
figure(1);
plot(xx,yy)
ylabel(txt{2})
xlabel(txt{1})

% Creating a range of velocities
v = linspace(0, 150, length(xx));

% Finding corresponding efficiency for every velocity configuration
N = eff( v, xx, yy, n, D );

% Plotting velocity vs efficiency
figure(2);
plot(v,N)
ylabel('Efficiency')
xlabel('Velocity')
