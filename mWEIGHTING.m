function [WEIGHTS] = mWEIGHTING(CATEGORIES, TYPE)
% Calculate a matrix of weights to be used while calculating agreement
%
%   CATEGORIES is a numeric vector specifying the possible categories.
%
%   TYPE is a parameter specifying the weighting scheme to be used for
%   partial agreement. The three options are below:
%       'identity' is for unordered/nominal categories
%       'linear' is for ordered categories and is relatively strict
%       'quadratic' is for ordered categories and is relatively forgiving
%
%   WEIGHTS is a q-by-q numeric matrix where q is the number of elements in
%   the CATEGORIES vector. Each value in this matrix specifies the amount
%   of "partial credit" that should be awarded for one rater assigning an
%   object to category k and another rater assigning that object to
%   category l, where k and l index rows and columns of WEIGHTS.
%
%   Example usage: mWEIGHTING(1:5, 'linear')
%   
%   (c) Jeffrey M Girard, 2016-2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate weights based on data scale
q = length(CATEGORIES);
if isnumeric(CATEGORIES)
    cat = CATEGORIES;
else
    cat = 1:q;
end
WEIGHTS = nan(q);
catmin = min(cat);
catmax = max(cat);
for k = 1:q
    for l = 1:q
        switch TYPE
            case 'identity'
                WEIGHTS = eye(q);
            case 'linear'
                if k==l
                    WEIGHTS(k,l) = 1;
                else
                    dist = abs(cat(k) - cat(l));
                    maxdist = catmax - catmin;
                    WEIGHTS(k,l) = 1 - (dist / maxdist);
                end
            case 'quadratic'
                if k==l
                    WEIGHTS(k,l) = 1;
                else
                    distsq = (cat(k) - cat(l)) ^ 2;
                    maxdistsq = (catmax - catmin) ^ 2;
                    WEIGHTS(k,l) = 1 - (distsq / maxdistsq);
                end
            otherwise
                error('Weighting must be identity, linear, or quadratic.');
        end
    end
end

end
