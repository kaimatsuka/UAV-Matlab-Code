function y = W_SC(W_TO)
% 
% DESCRIPTION: 
%   Nicolai surface controls weight estimate 
%
% INPUTS:
%   W_TO = weight at take-off (lb) 
%
% OUTPUTS:
%   W_SC = weight of control surfaces (lb)
%
% SOURCE:
%   Nicoli weight estimate equations (MAE 154A resource)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = 1.08*(W_TO^0.7);

end