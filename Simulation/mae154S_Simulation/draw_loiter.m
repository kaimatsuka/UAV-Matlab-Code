function [loiter_traj, heading] = draw_loiter(start_N,start_E,tar_N,tar_E,next_N,next_E,N)
% 
% DESCRIPTION:
%   This function draws line from start to near target, and once it
%   approaches target, draws N circles around the target to loiter.
%

    % design parameter
    d_space  = 200;  % (ft) maximum limit of spacing btw each point
    r_loiter = 1000; % (ft) loiter radius
    
    % approach to near target (call intermediate point
    heading_tar = atan2(tar_E-start_E,tar_N-start_N);
    d_int = sqrt((tar_N-start_N)^2+(tar_E-start_E)^2)-r_loiter;
    int_N = d_int*cos(heading_tar)+start_N; 
    int_E = d_int*sin(heading_tar)+start_E; 
    
    int_traj = draw_line(start_N,start_E,int_N,int_E);

    % loiter in circle for N times
    theta_start = pi + heading_tar; % theta = 0 at north, theta = 90 at east
    heading_next = atan2(next_E-tar_E,next_N-tar_N); % heading from target to next point (range: -pi to pi)
    diff = heading_next-theta_start; 
    if diff < 0
        diff = diff+2*pi;
    end
    theta_end   = theta_start+2*pi*N+diff;
    theta = [theta_start:d_space/r_loiter:theta_end]'; % Mx1 vector
    circle_traj = [cos(theta) sin(theta)]; % Mx2 vector
    circle_traj = circle_traj*r_loiter+ones(size(circle_traj,1),1)*[tar_N tar_E];
    
%     loiter_traj = int_traj;
%     heading = 0;

        loiter_traj = [int_traj;
                       circle_traj];
                   
    heading = atan2(loiter_traj(end,2)-loiter_traj(end-1,2),...
                    loiter_traj(end,1)-loiter_traj(end-1,1));
    
end