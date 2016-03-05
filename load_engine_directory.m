% engine_directory.m
%
% DESCRIPTION:
%   This directory stores all the information from all the engines we
%   consider
%
% INPUTS:
%   N/A
% 
% OUTPUTS:
%   This function is called by the optimization code to select the best
%   engine. 
%
% REVISION HISTORY:
%   02/20 File Created. 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eng_data = xlsread('Engine Directory.xlsx','EngInfo','B4:L12');
eng_names = cellstr(['RCG 15cc  ';'RCG 20cc  ';'RCG 26cc  ';'RCG 30cc  ';...
    'RCG 30cc_t';'RCG 55cc  ';'RCG 61cc  ';'NGH GT35  ';'NGH GT35R ']);

for i = 1:9
    engines(i).name     = eng_names(i); % engine names
    engines(i).P_max    = eng_data(i,1); % engine max power (hp)
    engines(i).eff      = eng_data(i,2); % engine efficiency
    engines(i).P_avail  = eng_data(i,3); % engine available power (hp)
    engines(i).d        = eng_data(i,4); % engine cylinder diameter (ft)
    engines(i).l        = eng_data(i,5); % engine cylinder length (ft)
    engines(i).W        = eng_data(i,6); % engine box width (ft)
    engines(i).L        = eng_data(i,7); % engine box length (ft)
    engines(i).H        = eng_data(i,8); % engine box depth (ft)
    engines(i).vol      = eng_data(i,9); % engine total volume (ft^3)
    engines(i).weight   = eng_data(i,10); % engine weight (lbs)
    engines(i).rpm      = eng_data(i,10); % engine rpm 
    
end

