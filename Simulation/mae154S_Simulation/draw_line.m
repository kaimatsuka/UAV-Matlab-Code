function [line_traj, heading] = draw_line(start_N,start_E,end_N,end_E)
% OUTPUT
%   line_traj = Mx2 matrix of trajectory
%   heading   = terminal heading direction of trajectory (rad: ranges -pi to pi)
%
    
    d_space = 200; % maximum limit of spacing btw each point
    d_total = sqrt((end_N-start_N)^2+(end_E-start_E)^2);
    N_p = ceil(d_total/d_space); % number of points
    line_traj = [linspace(start_N,end_N,N_p)' linspace(start_E,end_E,N_p)'];
    heading = atan2((end_N-start_N),(end_E-start_E));

end