% load_propeller.m
%
% DESCRIPTION:
%   This file load propeller data.
%
% INPUT:
%   none
%
% OUTPUT:
%
% REVISION HISTORY:
%   03/06 This file was created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [num,txt] = xlsread('Propeller Designnn.xlsx','Efficiency', 'A1:B35','basic');
% scale_factor_J     = 1
% scale_factor_eta_p = 1;

% [num,txt] = xlsread('Propeller Designnn.xlsx','Efficiency', 'C1:D44','basic');
% scale_factor_J     = 1
% scale_factor_eta_p = 1;

[num,txt] = xlsread('Propeller Designnn.xlsx','Efficiency', 'E1:F35','basic');
scale_factor_J     = 0.76; %scaling factor for small UAV
scale_factor_eta_p = 0.92; %scaling to take account of pusher prop

% Separating excel data, performing pchip interpolation

propeller.J     = num(:, 1)*scale_factor_J;
propeller.eta_p = num(:, 2)*scale_factor_eta_p;
