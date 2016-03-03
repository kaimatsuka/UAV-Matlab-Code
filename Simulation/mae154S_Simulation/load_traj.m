function [] = load_traj()

global traj_total

opt = 4;
switch opt 
    case 1 %---------------------------------------------------------------
        traj0 = [-300 0]; % starting point

        N1 = 1000;
        traj1 = [zeros(N1,1) linspace(0,35000,N1)'];
        traj1 = ones(size(traj1,1),1)*traj0(end,:)+traj1;
        N2 = 1000;
        traj2 = [zeros(N2,1) linspace(0,35000,N2)'];
        traj2 = ones(size(traj2,1),1)*traj1(end,:)+traj2;
        R = 1500; N3 = 1000;
        theta = linspace(0,pi,N3)';
        traj3 = [R*cos(theta)-R R*sin(theta)];
        traj3 = ones(size(traj3,1),1)*traj2(end,:)+traj3;
        N4 = 1000;
        traj4 = [zeros(N4,1) linspace(0,-35000,N4)'];
        traj4 = ones(size(traj4,1),1)*traj3(end,:)+traj4;

        traj_total = [traj1;
                      traj2
                      traj3;
                      traj4;];
        traj_total = fliplr(traj_total);
        
    case 2 %---------------------------------------------------------------
        % source: 
        %    http://stackoverflow.com/questions/13894715/draw-equidistant-points-on-a-spiral
        
        radius_start = 5000; % (ft)
        radius_end   = 35000; 
        thetaMax = 10*2*pi;
        awayStep = (radius_end-radius_start)/thetaMax;
        N = 1000;
        dist = linspace(0,2.2*10^6,N)';
        r    = linspace(radius_start, radius_end,N)';
        theta = dist./r;
%         r = awayStep*theta + radius_start;
        x = cos(theta).*r;
        y = sin(theta).*r;
        traj_total = [x y];
        traj_total = fliplr(traj_total);
        
    case 3 %---------------------------------------------------------------
%         traj0 = [-300 0]; % starting point
% 
%         N1 = 1000;
%         traj1 = [linspace(0,-17000,N1)' linspace(0,50000,N1)'];
%         traj1 = ones(size(traj1,1),1)*traj0(end,:)+traj1;
            
        [traj1 heading] = draw_line(0,-300,50000,-17000);
        % create spiral
        N = 10; % turns
        radius_end = 20000;
        radius_start = 1000;
        thetaMax = N*2*pi;
        chord = 1000;
        theta = 0;
        E_c = -radius_start;
        N_c = 0;

        max_num = 5000;
        pos = zeros(max_num,2);
        for ii = 1:max_num 
            awayStep = radius_end/thetaMax;
            away = awayStep*theta+radius_start;
            pos(ii,:) = [sin(theta) cos(theta)].*away+[N_c E_c];

            if theta >= thetaMax
                break;
            end
            theta = theta + chord/away;
        end
        traj2 = pos(1:ii,:);
        traj2 = ones(size(traj2,1),1)*traj1(end,:)+traj2;

        traj_total = [traj1;
                      traj2];
                      
    case 4 %---------------------------------------------------------------
        
        [traj1 heading] = draw_line(0,-300,50000,-17000);
        traj2 = draw_spiral(traj1(end,1),traj1(end,2), 10);
        
        [traj3 heading] = draw_line(traj2(end,1),traj2(end,2),55000,-10000);

        traj_total = [traj1;
                      traj2
                      traj3];
        
end


        

function [line_traj heading] = draw_line(start_N,start_E,end_N,end_E)
% OUTPUT
%   line_traj = Mx2 matrix of trajectory
%   heading   = terminal heading direction of trajectory (rad: ranges -pi to pi)
%
    
    d_space = 200; % maximum limit of spacing btw each point
    d_total = sqrt((end_N-start_N)^2+(end_E-start_E)^2);
    N_p = ceil(d_total/d_space); % number of points
    line_traj = [linspace(start_N,end_N,N_p)' linspace(start_E,end_E,N_p)'];
    heading = atan2((end_N-start_N),(end_E-start_E));


function spiral_traj = draw_spiral(start_N,start_E, N, psi)

    
    % Design parameters
    r_end   = 20000; % outer radius (ft)
    r_start = 1000;  % inner radius (ft)
    d_space = 200;   % spacing btw each points (ft)
    max_allowed_data = 10000; 
    
%     theta_start = -psi; 
%     theta_end = N*2*pi-psi; %(rad)
    theta_start = 0; 
    theta_end = N*2*pi; %(rad)
    E_c = start_E-r_start;
    N_c = start_N;

    %Initialization
    theta = theta_start;
    temp_pos = zeros(max_allowed_data,2);
    
    for ii = 1:max_allowed_data 
        awayStep = (r_end-r_start)/theta_end;
        away = awayStep*theta+r_start;
        temp_pos(ii,:) = [sin(theta) cos(theta)].*away+[N_c E_c];

        if theta >= theta_end
            break;
        end
        theta = theta + d_space/away;
    end
    
    spiral_traj = temp_pos(1:ii,:);