% Test Trajectory
clear all
close all
global traj_total

load_traj;

figure(2)
plot(traj_total(:,2),traj_total(:,1)'.')
hold on, grid on
axis('equal')