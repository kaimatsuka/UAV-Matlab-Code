function C_Dp_total = C_Dpi(K,Q,C_f,S_wet,S_ref,C_Dmisc,C_DLP)
% 
% DESCRIPTION: 
%   This function computes parasite drag 
%
% INPUTS:
%   K       = form factor (NxM) 
%   Q       = Intereference factor (Nx1)
%   C_f     = Intereference factor (NxM)
%   S_wet   = Intereference factor (Nx1)
%   S_ref   = Intereference factor (1x1)
%   C_Dmisc = miscellaneous drag coefficient (1x1)
%   C_DLP   = leakage drag coefficent (1x1)
%
% OUTPUTS:
%   C_Dp_total = parasite drag coefficient (1xM)
%
% where 
%   N = number of components, and
%   M = number of data points in velocity.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(K,1);
M = size(K,2);

C_Dpi = zeros(N,M); 

for i=1:N
    % for each component, calculate C_Dpi
    C_Dpi(i,:) = (K(i,:)*Q(i).*C_f(i,:)*S_wet(i))/S_ref;
    
end

% Sum drag components along 1st dimension, and add miscelaneous components
C_Dp_total = sum(C_Dpi,1)+ones(1,M)*(C_Dmisc+C_DLP); 

end