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


good_data  = data(idx1);
bad_data   = data(idx2);

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

    order_of_mag = order((max_val-min_val)); % find relevant order of magnitude
    max_closest = 10^order_of_mag*ceil(max_val/10^order_of_mag);  
    min_closest = 10^order_of_mag*floor(min_val/10^order_of_mag);

    % define edges
    edges = [min_closest:10^(order_of_mag):max_closest];

    N_good = histc(good_data, edges);
    N_good = fliplr(N_good); % for some reason, histc
    N_bad = histc(bad_data, edges);
    N_bad = fliplr(N_bad);

    N_total = [N_good; N_bad];

    axis_name = linspace(max_closest,min_closest,length(N_good)); % fix last arg
    axis_name = axis_name+(axis_name(1)-axis_name(2))/2;
    barh(axis_name,N_total',1,'stacked');
end
