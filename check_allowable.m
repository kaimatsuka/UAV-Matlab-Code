function var = check_allowable(original, sd)
% DESCRIPTION:
%   This file checks if any deviated lengths are negative (not possible)
%
% INPUT:
%   original: original variable
%   sd: how much deviation
%
% OUTPUT:
%   var: correct possible variable
%
% REVISION HISTORY:
%   02/27: File created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (original+sd) <= 0
        var = original;
    else
        var = original + sd;
    end
end