function [ALPHAK] = FAST_ALPHAK(CODES)
% Calculate Krippendorff's Alpha for dichotomous coding from two coders
%   [ALPHAK] = FAST_ALPHAK(CODES)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., coder).
%   This function can only handle two coders and two categories (use the
%   FULL_ALPHAK function for more than two coders and many categories).
%
%   ALPHAK is a chance-adjusted index of agreement. It estimates chance
%   agreement using a distribution-based approach. ALPHAK ranges from
%   -1 to 1 where 0 means coders were no better than chance.
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
[n,j] = size(CODES);
x = unique(CODES);
%% Output basic descriptives
fprintf('Number of items = %d\n',n);
fprintf('Number of coders = %d\n',j);
fprintf('Number of possible categories = %d\n',2);
fprintf('Observed categories = %s\n',mat2str(x));
%% Check for valid dichotomous data from two coders
if n < 1
    ALPHAK = NaN;
    fprintf('ERROR: At least 1 item is required.\n')
    return;
end
if j ~= 2
    ALPHAK = NaN;
    fprintf('ERROR: Use FULL_ALPHAK for more than 2 coders.\n');
    return;
end
if length(x) ~= 2
    ALPHAK = NaN;
    fprintf('ERROR: Use FULL_ALPHAK for more than 2 categories.\n');
    return;
end
%% Calculate contingency table and derivatives
a = sum(CODES(:,1)==x(2) & CODES(:,2)==x(2));
b = sum(CODES(:,1)==x(1) & CODES(:,2)==x(2));
c = sum(CODES(:,1)==x(2) & CODES(:,2)==x(1));
d = sum(CODES(:,1)==x(1) & CODES(:,2)==x(1));
m1 = (2*a + b + c) / 2;
m2 = (b + c + 2*d) / 2;
%% Calculate reliability and its components
P_O = (a + d) / n;
P_C = ((2 * m1) / (2 * n)) * ((2 * m1 - 1) / (2 * n - 1)) + ((2 * m2) / (2 * n)) * ((2 * m2 - 1) / (2 * n - 1));
ALPHAK = (P_O - P_C) / (1 - P_C);
%% Output reliability and its components
fprintf('Percent observed agreement = %.3f\n',P_O);
fprintf('Percent chance agreement = %.3f\n',P_C);
fprintf('\nKrippendorff''s Alpha = %.3f\n',ALPHAK);

end