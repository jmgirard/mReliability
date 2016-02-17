function [AC] = FAST_AC(CODES)
% Quickly calculate Gwet's AC1 for dichotomous coding from two coders
%   [AC] = FAST_AC(CODES)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., coder).
%   This function can only handle two coders and two categories (use the
%   FULL_AC function for more than two coders and many categories).
%
%   AC is a chance-adjusted index of agreement. In the dichotomous case, it
%   estimates chance agreement using a distribution-based approach. AC
%   ranges from -1 to 1 where 0 means coders were no better than chance.
%
%   Example usage: [AC] = FAST_AC(smiledata);
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
    AC = NaN;
    fprintf('ERROR: At least 1 item is required.\n')
    return;
end
if j ~= 2
    AC = NaN;
    fprintf('ERROR: Use FULL_AC for more than 2 coders.\n');
    return;
end
if length(x) ~= 2
    AC = NaN;
    fprintf('ERROR: Use FULL_AC for more than 2 categories.\n');
    return;
end
%% Calculate contingency table cells and marginals
a = sum(CODES(:,1)==x(2) & CODES(:,2)==x(2));
b = sum(CODES(:,1)==x(1) & CODES(:,2)==x(2));
c = sum(CODES(:,1)==x(2) & CODES(:,2)==x(1));
d = sum(CODES(:,1)==x(1) & CODES(:,2)==x(1));
m1 = (2*a + b + c) / 2;
m2 = (b + c + 2*d) / 2;
%% Calculate reliability and components
P_O = (a + d) / n;
P_C = (m1 / n) * (m2 / n) + (m1 / n) * (m2 / n);
AC = (P_O - P_C) / (1 - P_C);
%% Output reliability and variance components
fprintf('Percent observed agreement = %.3f\n',P_O);
fprintf('Percent chance agreement = %.3f\n',P_C);
fprintf('\nGwet''s AC = %.3f\n',AC);

end