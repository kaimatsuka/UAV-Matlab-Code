% load_airfoils.m
%
% DESCRIPTION:
%   This file load libraries of airfoil data
%
% INPUT:
%   none
%
% OUTPUT:
%   airfoils(ii).name
%               .alpha
%               .Cl
%               .Cd
%               .Cdp
%               .CM
%               .ClCd
%               .maxClCd
%               .Cl_alpha_deg
%               .Cl_alpha_rad
%               .alpha0
%
% REVISION HISTORY:
%   02/18: File created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%TODO: add library of maximum thickness ratio, chordwise location of the
%airfoil max thickness location, and chord-circumference ratio

airfoil_options = cellstr(['0008'; '0024'; '1408'; '1410'; '1412'; ...
                           '2412'; '2415'; '2418'; '4412'; '4415'; ...
                           '4418'; '6412']);
file_in = './Airfoil_Raw/';
airfoil_data = strcat(file_in, 'naca', airfoil_options, '.txt');
clrstr = cellstr(['-m';'-c';'-r';'-g';'-b';'-k';...
                ':m';':c';':r';':g';':b';':k']);

for ii = 1:length(airfoil_data)

    raw = dlmread(char(airfoil_data(ii)));
    
    airfoils(ii).name  = strcat('naca',airfoil_options(ii));
    airfoils(ii).alpha = raw(:,1);
    airfoils(ii).Cl    = raw(:,2);
    airfoils(ii).Cd    = raw(:,3);
    airfoils(ii).Cdp   = raw(:,4);
    airfoils(ii).CM    = raw(:,5);
    airfoils(ii).ClCd  = airfoils(ii).Cl./airfoils(ii).Cd;
    airfoils(ii).maxClCd = max(airfoils(ii).ClCd);
    
    % calculate max Cl, max alpha
    [maxCl maxInd] = max(airfoils(ii).Cl);
    airfoils(ii).maxCl    = maxCl;
    airfolis(ii).maxAlpha = airfoils(ii).alpha(maxInd);
    
    start_ind = find(airfoils(ii).alpha == -5);
    end_ind = find(airfoils(ii).alpha == 5.1);
    airfoils(ii).Cl_alpha_deg = (airfoils(ii).Cl(end_ind)-airfoils(ii).Cl(start_ind))/...
                            (airfoils(ii).alpha(end_ind)-airfoils(ii).alpha(start_ind));
    airfoils(ii).Cl_alpha_rad = airfoils(ii).Cl_alpha_deg*(180/pi);
    
    % Calculate the zero lift angle of attack --> alpha_0 (degrees)
    tmp = abs(airfoils(ii).Cl-0);                % find when Cl = 0
    [val idx] = min(tmp);                       % locate the index
    airfoils(ii).alpha0 = airfoils(ii).alpha(idx);% in degrees
    
end   