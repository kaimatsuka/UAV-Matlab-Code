function [spiral_traj, heading] = draw_spiral(start_N,start_E, psi)

    % Design parameters
    r_end   = 21000; % outer radius (ft) set 21000 to get 50 mile^2 coverage
    r_start = 1000;  % inner radius (ft)
    d_space = 200;   % spacing btw each points (ft)
    width   = 400;   % width between each spiral lines
    max_allowed_data = 10000; 
    
%     theta_start = -psi; 
%     theta_end = N*2*pi-psi; %(rad)
    N = (r_end-r_start)/d_space/20;
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
    heading = atan2(temp_pos(ii,2)-temp_pos(ii-1,2),temp_pos(ii,1)-temp_pos(ii-1,1));
    
end