function [ N ] = eff( v, xx, yy, n, D )

diff = zeros(length(v), length(xx));

% Nested for loop to find the difference of velocity to velocity obtained
% from parameter. Nested for loop ensures all combinations are considered
for j = 1:length(v) 
    
    for i = 1:length(xx)
        
        diff(j, i) = abs(v(j) - n*D*xx(i));
        
    end
end

% Sizing the matrix
[row, co] = size(diff);

N = zeros(length(yy));

% Finding the smallest difference in every column and saving this value.
% Then, finding index where this value occurs and finding the corresponding
% efficiency
for q = 1:row
    
        [Y, I] = min(diff(q, :));
        
        N(q) = yy(I);
        
end


end

