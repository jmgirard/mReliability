function [SA] = FULL_SA(CODES, CATEGORIES)
% Calculate specific agreement for each category
%   [SA] = FULL_SA(CODES, CATEGORIES)
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
%   SA is a vector containing specific agreement scores for each category.
%
%   Example usage: [SA] = FULL_SA(smiledata,[0,1]);
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

%% Remove items with all missing codes
CODES(all(~isfinite(CODES),2),:) = [];
%% Calculate basic descriptives
[n,r] = size(CODES);
x = unique(CODES);
x(~isfinite(x)) = [];
if nargin < 2
    CATEGORIES = x;
end
if exist('CATEGORIES','var')==0
    CATEGORIES = x;
end
if isempty(CATEGORIES)
    CATEGORIES = x;
end
CATEGORIES = unique(CATEGORIES(:));
q = length(CATEGORIES);
%% Output basic descriptives
fprintf('Number of items = %d\n',n);
fprintf('Number of raters = %d\n',r);
fprintf('Possible categories = %s\n',mat2str(CATEGORIES));
fprintf('Observed categories = %s\n\n',mat2str(x));
%% Check for valid data from more than one rater
if n < 1
    SA = NaN;
    fprintf('ERROR: At least 1 item is required.\n')
    return;
end
if r < 2
    SA = NaN;
    fprintf('ERROR: At least 2 raters are required.\n');
    return;
end
if any(ismember(x,CATEGORIES)==0)
    SA = NaN;
    fprintf('ERROR: Categories were observed in CODES that were not included in CATEGORIES.\n');
    return;
end
%% Calculate specific agreement for each category
SA = nan(1,q);
for k = 1:q
    numerator = 0;
    denominator = 0;
    for i = 1:n
        r_i = sum(isfinite(CODES(i,:)));
        if r_i < 2, continue; end
        r_ik = sum(CODES(i,:)==CATEGORIES(k));
        numerator = numerator + (r_ik * (r_ik - 1));
        denominator = denominator + (r_ik * (r_i - 1));
    end
    SA(k) = numerator / denominator;
    fprintf('Specific agreement for category %d = %.3f\n',CATEGORIES(k),SA(k));
end

end