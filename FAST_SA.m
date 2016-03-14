function [SA] = FAST_SA(CODES)
% Calculate specific agreement using simplified formulas
%   [SA] = FAST_SA(CODES)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., rater).
%   This function can handle two raters and two (dichotomous) categories.
%
%   SA is a vector containing specific agreement scores for each category.
%
%   Example usage: [SA] = FAST_SA(smiledata);
%   
%   (c) Jeffrey M Girard, 2016
%   
%   References:
%
%   Cicchetti, D. V., & Feinstein, A. R., High agreement but low kappa: II.
%   Resolving the paradoxes. Journal of Clinical Epidemiology, 1990, 43,
%   551-558.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remove items with any missing codes
CODES(any(~isfinite(CODES),2),:) = [];
%% Calculate basic descriptives
[n,r] = size(CODES);
x = unique(CODES);
x(~isfinite(x)) = [];
q = length(x);
%% Output basic descriptives
fprintf('Number of items = %d\n',n);
fprintf('Number of raters = %d\n',r);
fprintf('Observed categories = %s\n\n',mat2str(x));
%% Check for valid data from more than one rater
if n < 1
    SA = NaN;
    fprintf('ERROR: At least 1 item is required.\n')
    return;
end
if r ~= 2
    SA = NaN;
    fprintf('ERROR: Exactly 2 raters are required.\n');
    return;
end
if q > 2
    SA = NaN;
    fprintf('ERROR: Two or fewer categories are required.\n');
    return;
end
%% Calculate specific agreement for each category
n_11 = sum(CODES(:,1)==x(1) & CODES(:,2)==x(1));
n_22 = sum(CODES(:,1)==x(2) & CODES(:,2)==x(2));
n_12 = sum(CODES(:,1)==x(1) & CODES(:,2)==x(2));
n_21 = sum(CODES(:,1)==x(2) & CODES(:,2)==x(1));
SA(1) = (2 * n_11) / (2 * n_11 + n_12 + n_21);
SA(2) = (2 * n_22) / (2 * n_22 + n_12 + n_21);
fprintf('Specific agreement for category %d = %.3f\n',x(1),SA(1));
fprintf('Specific agreement for category %d = %.3f\n',x(2),SA(2));

end