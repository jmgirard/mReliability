function [SA] = FULL_SA(CODES, CATEGORIES)
%   (c) Jeffrey M Girard, 2015
%   
%   Reference: http://www.john-uebersax.com/stat/raw.htm

%% Remove items with all missing codes
CODES(all(~isfinite(CODES),2),:) = [];
%% Calculate basic descriptives
[n,r] = size(CODES);
x = unique(CODES);
x(~isfinite(x)) = [];
if nargin < 2
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
        r_ik = sum(CODES(i,:)==CATEGORIES(k));
        numerator = numerator + (r_ik * (r_ik - 1));
        denominator = denominator + (r_ik * (r_i - 1));
    end
    SA(k) = numerator / denominator;
    fprintf('Specific agreement for category %d = %.3f\n',CATEGORIES(k),SA(k));
end

end