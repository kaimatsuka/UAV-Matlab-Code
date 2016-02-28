function CL = velocity2CL(velocity,weight,alt,n)
% velocity2CL.m
%
% DESCRIPTION:
%   This function calculates CL for given velocity.
%
% INPUT:
%   velocity = velocity of interest (ft/s)
%   weight   = weight of UAV at this time (lb)
%   alt      = altitude of interest (ft)
%   n        = load factor (N/A)
%
% OUTPUT:
%
% REVISION HISTORY:
%   02/28: First file is created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    error('Error: Not enough inputs')
elseif nargin < 4
    n = 1; % default value
end
    
rho = calc_atmos(alt);
CL = weight*n/(0.5*rho*velocity^2);

end