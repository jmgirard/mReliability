function [ALPHAK, P_O, P_C] = mALPHAK(CODES, CATEGORIES, WEIGHTING)
% Calculate Krippendorff's alpha coefficient using generalized formulas
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
%   ALPHAK is a chance-adjusted index of agreement.
%
%   P_O is the percent observed agreement (from 0.000 to 1.000).
%
%   P_C is the estimated percent chance agreement (from 0.000 to 1.000).
%
%   Example usage: mALPHAK(fishdata, [1, 2, 3], 'identity')
%   
%   (c) Jeffrey M Girard, 2016-2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remove items that do not have codes from at least two raters
CODES(sum(isfinite(CODES),2)<2,:) = [];
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
    ALPHAK = NaN;
    fprintf('ERROR: At least 1 valid object is required. \n')
    return;
end
if r < 2
    ALPHAK = NaN;
    fprintf('ERROR: At least 2 raters are required. \n');
    return;
end
if any(ismember(x, CATEGORIES) == 0)
    ALPHAK = NaN;
    fprintf('ERROR: Unexpected category in CODES. \n');
    return;
end
%% Get weights from mWEIGHTING function
weights = mWEIGHTING(CATEGORIES, WEIGHTING);
%% Create n-by-q matrix (rater counts in item by category matrix)
r_ik = zeros(n, q);
for k = 1:q
    codes_k = CODES == CATEGORIES(k);
    r_ik(:,k) = codes_k * ones(r,1);
end
rstar_ik = transpose(weights * transpose(r_ik));
%% Calculate percent observed agreement
r_i = r_ik * ones(q, 1);
r_ik = r_ik(r_i >= 2, :);
rstar_ik = rstar_ik(r_i >= 2, :);
r_i = r_i(r_i >= 2);
rbar_i = mean(r_i);
nprime = size(r_ik, 1);
epsilon = 1 / sum(r_i);
observed = (r_ik .* (rstar_ik - 1)) * ones(q, 1);
possible = rbar_i .* (r_i - 1);
P_O = (1 - epsilon) .* sum(observed ./ (possible)) ./ nprime + epsilon;
%% Calculate percent chance agreement
pihat = transpose(repmat(1 / n, 1, n) * (r_ik ./ rbar_i));
P_C = sum(sum(weights .* (pihat * transpose(pihat))));
%% Calculate reliability point estimate
ALPHAK = (P_O - P_C) / (1 - P_C);

end
