% airfoil_analysis.m
%
% DESCRIPTION:
%   This function analyzes raw airfoil data from Xfoil and converts the 2D
%   data into 3D lift coefficients
% 
% OUTPUTS:
%   This function outputs a structure with all the airfoils analyzed.
%
% REVISION HISTORY:
%   02/12 v 1.0
%   02/17 v 2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

airfoil_txt = cellstr(['0008'; '0024'; '1408'; '1410'; '1412'; ...
    '2412'; '2415'; '2418'; '4412'; '4415'; '4418'; '6412']);
file_in = './Airfoil_Raw/';
airfoil_data = strcat(file_in, 'naca', airfoil_txt, '.txt');
clrstr = cellstr(['-m';'-c';'-r';'-g';'-b';'-k';...
                ':m';':c';':r';':g';':b';':k']);
            
% Airfoil(ii) is airfoil listed in airfoils vector
% Text file is formatted in columns:
% Alpha (deg) CL CD CDp CM Top_Xtr Bot_Xtr
for ii = 1:length(airfoil_data)
% for ii = 1:1
    raw = dlmread(char(airfoil_data(ii)));
    
    airfoil(ii).name  = strcat('naca',airfoil_txt(ii));
    airfoil(ii).alpha = raw(:,1);
    airfoil(ii).Cl    = raw(:,2);
    airfoil(ii).Cd    = raw(:,3);
    airfoil(ii).Cdp   = raw(:,4);
    airfoil(ii).CM    = raw(:,5);
    airfoil(ii).ClCd  = airfoil(ii).Cl./airfoil(ii).Cd;
    
    % Determine minimum drag condition; Cl/Cd max
    airfoil(ii).maxClCd = max(airfoil(ii).ClCd);
    
    % Plot Cl vs. Alpha graph (airfoil)
    figure(1);
    hold on; grid on;
    plot(airfoil(ii).alpha,airfoil(ii).Cl,char(clrstr(ii)));
    
    % Plot Cd vs. Alpha graph (airfoil)
    figure(2);
    hold on; grid on;
    plot(airfoil(ii).alpha, airfoil(ii).Cd, char(clrstr(ii)));
    
    % Plot Cl vs. Cd graph (airfoil)
    figure(3);
    hold on; grid on;
    plot(airfoil(ii).Cd, airfoil(ii).Cl, char(clrstr(ii)));
    
    % Plot Cl/Cd vs. Alpha graph (airfoil)
    figure(4);
    hold on; grid on;
    plot(airfoil(ii).alpha, airfoil(ii).ClCd, char(clrstr(ii)));
    
    legendInfo{ii} = strcat('NACA', char(airfoil_txt(ii)));
end

figure(1);
legend(legendInfo,'Location','Southeast');
xlabel('\alpha (degrees)'); ylabel('C_l');
title('C_l vs. \alpha of Select Airfoils')

figure(2);
legend(legendInfo,'Location','Southeast');
xlabel('\alpha (degrees)'); ylabel('C_d');
xlim([-10 10]);  ylim([0 0.02]);
title('C_d vs. \alpha of Select Airfoils')

figure(3);
legend(legendInfo,'Location','Northwest');
xlabel('C_d'); ylabel('C_l');
xlim([0 0.02]);
title('C_l vs. C_d of Select Airfoils')

figure(4);
legend(legendInfo, 'Location', 'Northwest');
xlabel('\alpha (degrees)'); ylabel('C_l/C_d');
title('C_l/C_d vs. \alpha of Select Airfoils')

%% calculate Cl_alpha & CL_alpha
load_UAV_parameters
for ii = 1:length(airfoil_data)
% for ii = 1:1
    % Determine indices where AoA is -5 and 5 degrees; where all airfoil
    % have straight line data
    start_ind = find(airfoil(ii).alpha == -5);
    end_ind   = find(airfoil(ii).alpha == 5.1);
    
    % Calculate the slope both /deg and /rad --> Cl_alpha
    airfoil(ii).Cl_alpha_deg = (airfoil(ii).Cl(end_ind)-airfoil(ii).Cl(start_ind))/...
                            (airfoil(ii).alpha(end_ind)-airfoil(ii).alpha(start_ind));
    airfoil(ii).Cl_alpha_rad = airfoil(ii).Cl_alpha_deg*(180/pi);
    
    % Calculate the zero lift angle of attack --> alpha_0 (degrees)
    tmp = abs(airfoil(ii).Cl-0);                % find when Cl = 0
    [idx idx] = min(tmp);                       % locate the index
    airfoil(ii).alpha0 = airfoil(ii).alpha(idx);% in degrees
    
    %Calculate the 3D lift curve slope --> CL_alpha
    [airfoil(ii).CL_alpha_rad, airfoil(ii).CL_alpha_deg] = CL_alpha(airfoil(ii).Cl_alpha_rad, wing.A, wing.e);
    
    %Calculate the 3D CL values
    airfoil(ii).CL = airfoil(ii).CL_alpha_deg*(airfoil(ii).alpha-airfoil(ii).alpha0);

    figure(1);
    hold on; grid on;
    plot(airfoil(ii).alpha, airfoil(ii).CL, char(clrstr(ii)));
    
    % Find CL and alpha corresponding to max Cl/Cd value
    idx = find(airfoil(ii).ClCd - airfoil(ii).maxClCd == 0);
    airfoil(ii).CL_maxClCd = airfoil(ii).CL(idx);
    airfoil(ii).alpha_maxClCd = airfoil(ii).alpha(idx);
end