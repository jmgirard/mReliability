function [ALPHAK] = FAST_ALPHAK(CODES)
% Calculate Krippendorff's Alpha for dichotomous coding from two raters
%   [ALPHAK] = FAST_ALPHAK(CODES)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., rater).
%   This function can only handle two raters and two categories (use the
%   FULL_ALPHAK function for more than two raters and many categories).
%
%   ALPHAK is a chance-adjusted index of agreement. It estimates chance
%   agreement using a distribution-based approach. ALPHAK ranges from
%   -1 to 1 where 0 means raters were no better than chance.
%
%   Example usage: [ALPHAK] = FAST_ALPHAK(smiledata);
%
%   (c) Jeffrey M Girard, 2016
%   
%   References:
%
%   Krippendorff, K. (1980). Content analysis: An introduction to its
%   methodology (1st ed.). Newbury Park, CA: Sage Publications.
%
%   Hayes, A. F., & Krippendorff, K. (2007). Answering the call for a
%   standard reliability measure for coding data. Communication Methods
%   and Measures, 1(1), 77–89.
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
    ALPHAK = NaN;
    fprintf('ERROR: At least 1 item is required.\n')
    return;
end
if r ~= 2
    ALPHAK = NaN;
    fprintf('ERROR: Use FULL_ALPHAK for more than 2 raters.\n');
    return;
end
if length(x) ~= 2
    ALPHAK = NaN;
    fprintf('ERROR: Use FULL_ALPHAK for more than 2 categories.\n');
    return;
end
%% Calculate contingency table and derivatives
n_22 = sum(CODES(:,1)==x(2) & CODES(:,2)==x(2));
n_12 = sum(CODES(:,1)==x(1) & CODES(:,2)==x(2));
n_21 = sum(CODES(:,1)==x(2) & CODES(:,2)==x(1));
n_11 = sum(CODES(:,1)==x(1) & CODES(:,2)==x(1));
m1 = (2*n_22 + n_12 + n_21) / 2;
m2 = (n_12 + n_21 + 2*n_11) / 2;
%% Calculate reliability and its components
P_O = (n_22 + n_11) / n;
P_C = ((2 * m1) / (2 * n)) * ((2 * m1 - 1) / (2 * n - 1)) + ((2 * m2) / (2 * n)) * ((2 * m2 - 1) / (2 * n - 1));
ALPHAK = (P_O - P_C) / (1 - P_C);
%% Output reliability and its components
fprintf('Percent observed agreement = %.3f\n',P_O);
fprintf('Percent chance agreement = %.3f\n',P_C);
fprintf('\nKrippendorff''s Alpha = %.3f\n',ALPHAK);

end