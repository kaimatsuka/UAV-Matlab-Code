function y = W_ES(W_FS,W_TRON)
% 
% DESCRIPTION: 
%   Nicolai electrical system weight estimate 
%
% INPUTS:
%   W_FS   = fuel system weight (lb) 
%   W_TRON = avionics/electronics bare weight (lb)
%
% OUTPUTS:
%   W_ES = weight of electrical system (lb)
%
% SOURCE:
%   Nicoli weight estimate equations (MAE 154A resource)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = 426*(((W_FS+W_TRON)/10^3)^0.51);

end