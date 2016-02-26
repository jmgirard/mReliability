function [GAMMA] = FAST_GAMMA(CODES)
% Calculate Gwet's gamma coefficient using simplified formulas
%   [GAMMA] = FAST_AC(CODES)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., rater).
%   This function can only handle two raters and two categories (use the
%   FULL_AC function for more than two raters and many categories).
%
%   GAMMA is a chance-adjusted index of agreement.
%
%   Example usage: [GAMMA] = FAST_GAMMA(smiledata);
%
%   (c) Jeffrey M Girard, 2016
%   
%   References:
%
%   Gwet, K. L. (2008). Computing inter-rater reliability and its variance
%   in the presence of high agreement. The British Journal of Mathematical
%   and Statistical Psychology, 61(1), 29–48.
%
%   Gwet, K. L. (2014). Handbook of inter-rater reliability: The definitive
%   guide to measuring the extent of agreement among raters (4th ed.).
%   Gaithersburg, MD: Advanced Analytics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remove items with any missing codes
CODES(any(~isfinite(CODES),2),:) = [];
%% Calculate basic descriptives
[n,r] = size(CODES);
x = unique(CODES);
%% Output basic descriptives
fprintf('Number of items = %d\n',n);
fprintf('Number of raters = %d\n',r);
fprintf('Number of possible categories = %d\n',2);
fprintf('Observed categories = %s\n',mat2str(x));
%% Check for valid dichotomous data from two raters
if n < 1
    GAMMA = NaN;
    fprintf('ERROR: At least 1 item is required.\n')
    return;
end
if r ~= 2
    GAMMA = NaN;
    fprintf('ERROR: Use FULL_GAMMA for more than 2 raters.\n');
    return;
end
if length(x) ~= 2
    GAMMA = NaN;
    fprintf('ERROR: Use FULL_GAMMA for more than 2 categories.\n');
    return;
end
%% Calculate contingency table cells and marginals
n_22 = sum(CODES(:,1)==x(2) & CODES(:,2)==x(2));
n_12 = sum(CODES(:,1)==x(1) & CODES(:,2)==x(2));
n_21 = sum(CODES(:,1)==x(2) & CODES(:,2)==x(1));
n_11 = sum(CODES(:,1)==x(1) & CODES(:,2)==x(1));
m1 = (2*n_22 + n_12 + n_21) / 2;
m2 = (n_12 + n_21 + 2*n_11) / 2;
%% Calculate reliability and components
P_O = (n_22 + n_11) / n;
P_C = (m1 / n) * (m2 / n) + (m1 / n) * (m2 / n);
GAMMA = (P_O - P_C) / (1 - P_C);
%% Output reliability and variance components
fprintf('Percent observed agreement = %.3f\n',P_O);
fprintf('Percent chance agreement = %.3f\n',P_C);
fprintf('\nGwet''s gamma coefficient = %.3f\n',GAMMA);

end