
% define discrete points on trajectory 

%
clear all
close all
global traj_total

opt = 1;
switch opt 
    case 1

        load_traj;
        figure(2)
        plot(traj_total(:,2),traj_total(:,1)'.')
        hold on, grid on
        plot(traj_total(end,2),traj_total(end,1),'xr')
%         plot(path1(:,1),path1(:,2),'xr')
        axis('equal')
        
    case 2
        in = [1 0 0 0 120];
        out = my_guid(in);
        in = [out(1) 0 120 -30  120];
        out = my_guid(in);
end
% current position
%{
traj1 = [ones(100,1)*-1000 linspace(0,35000,100)'];
traj2 = [zeros(100,1) linspace(0,35000,100)'];
traj2 = ones(size(traj2,1),1)*traj1(end,:)+traj2;
R = 1500; N = 100;
theta = linspace(0,pi,N)';
traj3 = [R*cos(theta)-R R*sin(theta)];
traj3 = ones(size(traj3,1),1)*traj2(end,:)+traj3;


traj_total = [traj1; 
              traj2;
              traj3];


N = 100;
path1 = [linspace(0,1000,N)' linspace(0,30000,N)'];


%
L1 = 2000; % ft
          
% initiate          
ind = 1;

% simulation
for ii = 1:1
    pos = path1(ii,:);
    traj_remain = traj_total(ind:end,:);
    x = (traj_remain-ones(size(traj_remain,1),1)*pos);
    d = (x(:,1).^2 + x(:,2).^2).^0.5;
    ind_remain = find(d>L1);
    ind = ind_remain(1);
end


    
          

%}

