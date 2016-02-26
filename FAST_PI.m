function [PI] = FAST_PI(CODES)
% Quickly calculate Scott's Pi for dichotomous coding from two raters
%   [PI] = FAST_PI(CODES)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., rater).
%   This function can only handle two raters and two categories (use the
%   FULL_PI function for more than two raters and many categories).
%
%   PI is a chance-adjusted index of agreement. It estimates chance
%   agreement using a distribution-based approach. PI ranges from -1 to 1
%   where 0 means raters were no better than chance.
%
%   Example usage: [PI] = FAST_PI(smiledata);
%
%   (c) Jeffrey M Girard, 2016
%   
%   References:
%
%   Scott, W. A. (1955). Reliability of content analysis: The case of
%   nominal scaling. Public Opinion Quarterly, 19(3), 321–325.
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
    PI = NaN;
    fprintf('ERROR: At least 1 item is required.\n')
    return;
end
if r ~= 2
    PI = NaN;
    fprintf('ERROR: Use FULL_PI for more than 2 raters.\n');
    return;
end
if length(x) ~= 2
    PI = NaN;
    fprintf('ERROR: Use FULL_PI for more than 2 categories.\n');
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
P_C = (m1 / n) * (m1 / n) + (m2 / n) * (m2 / n);
PI = (P_O - P_C) / (1 - P_C);
%% Output reliability and its components
fprintf('Percent observed agreement = %.3f\n',P_O);
fprintf('Percent chance agreement = %.3f\n',P_C);
fprintf('\nScott''s Pi = %.3f\n',PI);

end