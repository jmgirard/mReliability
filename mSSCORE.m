function [SSCORE, P_O, P_C] = mSSCORE(CODES, CATEGORIES, WEIGHTING)
% Calculate Bennett et al.'s S score using generalized formula
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single object of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., rater).
%   This function can handle any number of raters and values.
%
%   CATEGORIES is an optional parameter specifying the possible categories
%   as a numerical vector. If this variable is not specified, then the
%   possible categories are inferred from the CODES matrix. This can
%   underestimate reliability if all possible categories aren't used.
%
%   WEIGHTING is an optional parameter specifying the weighting scheme to
%   be used for partial agreement. The three options are below:
%       'identity' is for unordered/nominal categories (default)
%       'linear' is for ordered categories and is relatively strict
%       'quadratic' is for ordered categories and is relatively forgiving
%
%   SSCORE is a chance-corrected index of agreement. It assumes that each
%   category has an equal chance of being selected at random. It ranges
%   from -1.0* to 1.0 where 0.0 means raters were no better than chance.
%   *The actual lower bound is determined by the number of categories.
%
%   P_O is the percent observed agreement (from 0.000 to 1.000).
%
%   P_C is the estimated percent chance agreement (from 0.000 to 1.000).
%
%   Example usage: mSSCORE(fishdata, [1, 2, 3], 'linear')
%   
%   (c) Jeffrey M Girard, 2016-2018
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remove items with all missing codes
CODES(all(~isfinite(CODES), 2), :) = [];
%% Calculate basic descriptives
[n, r] = size(CODES);
x = unique(CODES);
x(~isfinite(x)) = [];
if nargin < 2
    CATEGORIES = x;
    WEIGHTING = 'identity';
elseif nargin < 3
    WEIGHTING = 'identity';
end
if isempty(CATEGORIES)
    CATEGORIES = x;
end
CATEGORIES = unique(CATEGORIES(:));
q = length(CATEGORIES);
%% Check for valid data from multiple raters
if n < 1
    SSCORE = NaN;
    fprintf('ERROR: At least 1 valid object is required. \n')
    return;
end
if r < 2
    SSCORE = NaN;
    fprintf('ERROR: At least 2 raters are required. \n');
    return;
end
if any(ismember(x, CATEGORIES) == 0)
    SSCORE = NaN;
    fprintf('ERROR: Unexpected category in CODES. \n');
    return;
end
%% Get weights from mWEIGHTING function
weights = mWEIGHTING(CATEGORIES, WEIGHTING);
%% Create n-by-q matrix (rater counts in item by category matrix)
r_ik = zeros(n, q);
for k = 1:q
    codes_k = CODES == CATEGORIES(k);
    r_ik(:,k) = codes_k * ones(r, 1);
end
rstar_ik = transpose(weights * transpose(r_ik));
%% Calculate percent observed agreement
r_i = r_ik * ones(q, 1);
observed = (r_ik .* (rstar_ik - 1)) * ones(q,1);
nprime = sum(r_i >= 2);
possible = r_i .* (r_i - 1);
P_O = sum(observed(r_i >= 2) ./ (possible(r_i >= 2))) ./ nprime;
%% Calculate percent chance agreement
P_C = sum(weights(:)) / (q ^ 2);
%% Calculate reliability point estimate
SSCORE = (P_O - P_C) / (1 - P_C);

end