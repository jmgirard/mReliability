function [KAPPA] = FAST_KAPPA(CODES)
% Quickly calculate Cohen's Kappa for dichotomous coding from two raters
%   [KAPPA] = FAST_KAPPA(CODES)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., rater).
%   This function can only handle two raters and two categories (use the
%   FULL_PI function for more than two raters and many categories).
%
%   KAPPA is a chance-adjusted index of agreement. It estimates chance
%   agreement using a distribution-based approach. KAPPA ranges from
%   -1 to 1 where 0 means raters were no better than chance.
%
%   Example usage: [KAPPA] = FAST_KAPPA(smiledata);
%
%   (c) Jeffrey M Girard, 2016
%   
%   References:
%
%   Cohen, J. (1960). A coefficient of agreement for nominal scales.
%   Educational and Psychological Measurement, 20(1), 37–46.
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
    KAPPA = NaN;
    fprintf('ERROR: At least 1 item is required.\n')
    return;
end
if r ~= 2
    KAPPA = NaN;
    fprintf('ERROR: Use FULL_KAPPA for more than 2 raters.\n');
    return;
end
if length(x) ~= 2
    KAPPA = NaN;
    fprintf('ERROR: Use FULL_KAPPA for more than 2 categories.\n');
    return;
end
%% Calculate contingency table and derivatives
n_22 = sum(CODES(:,1)==x(2) & CODES(:,2)==x(2));
n_12 = sum(CODES(:,1)==x(1) & CODES(:,2)==x(2));
n_21 = sum(CODES(:,1)==x(2) & CODES(:,2)==x(1));
n_11 = sum(CODES(:,1)==x(1) & CODES(:,2)==x(1));
f1 = n_22 + n_21;
g1 = n_22 + n_12;
f2 = n_12 + n_11;
g2 = n_21 + n_11;
%% Calculate reliability and its components
P_O = (n_22 + n_11) / n;
P_C = (f1 / n) * (g1 / n) + (f2 / n) * (g2 / n);
KAPPA = (P_O - P_C) / (1 - P_C);
%% Output reliability and its components
fprintf('Percent observed agreement = %.3f\n',P_O);
fprintf('Percent chance agreement = %.3f\n',P_C);
fprintf('\nCohen''s Kappa = %.3f\n',KAPPA);

end