function [SA] = mSPECIFIC(CODES, CATEGORIES, WEIGHTING)
% Calculate specific agreement for each category using generalized formula
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., rater).
%   This function can handle multiple raters and multiple categories.
%
%   CATEGORIES is an optional parameter specifying the possible categories
%   as a numerical vector. If this variable is not specified, then the
%   possible categories are inferred from CODES (in numerical order).
%
%   WEIGHTING is an optional parameter specifying the weighting scheme to
%   be used for partial agreement. The three options are below:
%       'identity' is for unordered/nominal categories (default)
%       'linear' is for ordered categories and is relatively strict
%       'quadratic' is for ordered categories and is relatively forgiving
%
%   SA is a vector containing specific agreement scores for each category.
%
%   Example usage: mSPECIFIC(fishdata, [1, 2, 3], 'identity')
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
    SA = NaN;
    fprintf('\n ERROR: At least 1 valid object is required. \n')
    return;
end
if r < 2
    SA = NaN;
    fprintf('\n ERROR: At least 2 raters are required. \n');
    return;
end
if any(ismember(x, CATEGORIES) == 0)
    SA = NaN;
    fprintf('\n ERROR: Unexpected category in CODES. \n');
    return;
end
%% Get weights from mWEIGHTING function
weights = mWEIGHTING(CATEGORIES, WEIGHTING);
%% Create n-by-q matrix (rater counts in item by category matrix)
r_ik = zeros(n, q);
for k = 1:q
    codes_k = CODES == CATEGORIES(k);
    r_ik(:, k) = codes_k * ones(r, 1);
end
%% Weight n-by-q matrix based on weighting scheme
rstar_ik = transpose(weights * transpose(r_ik));
%% Calculate category-specific agreement
r_i = r_ik * ones(q,1);
numerator = sum(r_ik .* (rstar_ik - 1));
denominator = sum(r_ik .* (r_i - 1));
SA = numerator ./ denominator;
%% Output reliability and variance components
fprintf('Specific agreement for category %d = %.3f \n', [CATEGORIES'; SA]);

end