function [dist_diff] = latlon2NE(latlon_start,latlon_end,lat_ref)
% DESCRIPTION:
%   This function computes the relative distance between two latitude and
%   longitude to the distance difference in north east.
%
% INPUT:
%   latlon_start = latitude and longitude of starting point (deg) 1x2
%   latlon_end   = latitude and longitudeof starting point (deg)  1x2
%   lat_ref   = latitude of reference point (deg) This should be
%       approximately center of whole trajectory. This is used in 
%       caluclation of dist_E. 
%
% OUTPUT:
%   dist_N = distance in North direction from start to end
%   dist_E = distance in East direction from start to end
%

    % Earth Radius
    mile2ft = 5280;
    E_R = 3959*mile2ft; 
    
    dist_N = E_R*(latlon_end(1)-latlon_start(1))*pi/180;
    dist_E = E_R*cosd(lat_ref)*(latlon_end(2)-latlon_start(2))*pi/180;
    
    dist_diff = [dist_N dist_E];
end
