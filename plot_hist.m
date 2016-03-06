function [] = plot_hist(data,idx1,idx2)
% plot_hist.m
%
% DESCRIPTION:
%   This funciton plots a historgram from 'data' passed down. It stacks the
%   part of data specified by index1 and put data specified by index2 on
%   top.
%
% INPUT:
%   data  = data of interest  (1xM)
%   indx1 = an array of indices of good data (1xN)
%   indx2 = an array of indices of bad data  (1xM-N)
%
% OUTPUT:
%   plot of Histogram
%
% REVISION HISTORY:
%   02/23: file created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


max_val = max(data); %find max val of data (wiht as many sig figs)
min_val = min(data); %find max val of data (wiht as many sig figs)


if max_val == min_val
    % TODO: work on here
    warning('Histogram cannot be plotted because data is constant');
    %{
    N_good = length(good_data);
    N_bad  = length(bad_data);
    N_total = [N_good; N_bad];
    
    barh(good_data(1),N_total,1,'stacked');
    %}
else
    
    % find closest value of relevant sig figs.
    %
    %   ex. if max_val = 15.432 and min_val = 8.67, relevant order of magnitude
    %       is set to be 0.1. Then following will calculate max_closest = 15.5,
    %       and min_closest = 8.6
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N_edges = 20;
    
    PERCENT_SUCCESS = 0.2;
    
    order_of_mag = order((max_val-min_val))-1; % find relevant order of magnitude
    max_closest = 10^order_of_mag*ceil(max_val/10^order_of_mag);  
    min_closest = 10^order_of_mag*floor(min_val/10^order_of_mag);
    
    increment = (max_closest-min_closest)/N_edges;

    % define edges
    edges = [min_closest:increment:max_closest];
    
    if idx1(1) ~= -1
        %TODO: incorporate case where fuel is 
        pass_data  = data(idx1);
        N_pass = histc(pass_data, edges);
        N_pass = fliplr(N_pass); % for some reason, histc
        
%         [sorted_weight sort_ind] = sort(arrayfun(@(x) x.weight.total, pass_data)); 
%         ind_pass_good = sort_ind(1:NUM_ITERATION/10);     % select good aircraft (lowest 10% in weight)
%         ind_pass_bad  = sort_ind(NUM_ITERATION/10+1:end); % rest of aircraft
        
        
    else
        N_pass = zeros(length(edges));
    end

    if idx2(1) ~= -1
        %TODO: incorporate case where fuel is 
        fail_data   = data(idx2);
        N_fail = histc(fail_data, edges);
        N_fail = fliplr(N_fail);
    else
        N_fail = zeros(length(edges));
    end

    N_total = [N_pass; N_fail];

    axis_name = linspace(max_closest,min_closest,length(N_pass)); % fix last arg
    axis_name = axis_name+(axis_name(1)-axis_name(2))/2;
    barh(axis_name,N_total',1,'stacked');
end
