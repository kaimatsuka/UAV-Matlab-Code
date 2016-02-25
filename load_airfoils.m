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
%   02/20: Edit to include max thickness, location of max thickness,
%   perimeter, CL/Clmax values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

airfoil_options = cellstr(['0008'; '0024'; '1408'; '1410'; '1412'; ...
                           '2412'; '2415'; '2418'; '4412'; '4415'; ...
                           '4418'; '6412']);
file_in = './Airfoil_Raw/';
airfoil_data = strcat(file_in, 'naca', airfoil_options, '.txt');
airfoil_geo_dat = strcat(file_in, 'naca', airfoil_options, '_geo.txt');
max_thick = [0.300; 0.300; 0.300; 0.300; 0.299;...
    0.300; 0.295; 0.300; 0.300; 0.309; 0.300; 0.301]; 
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
    airfoils(ii).Clmax = max(airfoils(ii).Cl);
    idx = find(airfoils(ii).Clmax - airfoils(ii).Cl == 0);
    airfoils(ii).Clmax_alpha = airfoils(ii).alpha(idx);
    airfoils(ii).ClCd  = airfoils(ii).Cl./airfoils(ii).Cd;
    airfoils(ii).maxClCd = max(airfoils(ii).ClCd);

    %%%%%%%%%%%%%% GEOMETERIC AIRFOIL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE: All GEOMETRIC PROPERTIES will be in its own structure
    airfoils(ii).GEO.max_thick          = str2num(airfoil_options{ii}(3:4))/100;
    airfoils(ii).GEO.max_thick_location = max_thick(ii);
    
    %PERIMETER CALCULATIONS
    raw_geo = dlmread(char(airfoil_geo_dat(ii)));
    perimeter = 0;
    perimeter = perimeter + norm(raw_geo(1:end-1,:) - raw_geo(2:end,:));
    perimeter = perimeter + norm(raw_geo(end,:)- raw_geo(1,:));
    airfoils(ii).GEO.perimeter = perimeter;
    
    % Calculate the zero lift angle of attack --> alpha_0 (degrees)
    tmp = abs(airfoils(ii).Cl-0);                % find when Cl = 0
    [idx idx] = min(tmp);                       % locate the index
    airfoils(ii).alpha0 = airfoils(ii).alpha(idx);% in degrees
    
    % calculate max Cl, max alpha
    [maxCl maxInd] = max(airfoils(ii).Cl);
    airfoils(ii).maxCl    = maxCl;
    airfoils(ii).maxAlpha = airfoils(ii).alpha(maxInd);
    
    % Determine indices where AoA is -5 and 5 degrees; where all airfoil
    % have straight line data
    start_ind = find(airfoils(ii).alpha == -5);
    end_ind = find(airfoils(ii).alpha == 5.1);
    
    % Calculate the slope both /deg and /rad --> Cl_alpha
    airfoils(ii).Cl_alpha_deg = (airfoils(ii).Cl(end_ind)-airfoils(ii).Cl(start_ind))/...
                            (airfoils(ii).alpha(end_ind)-airfoils(ii).alpha(start_ind));
    airfoils(ii).Cl_alpha_rad = airfoils(ii).Cl_alpha_deg*(180/pi);
end   