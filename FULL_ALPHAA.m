function [ ALPHAA ] = FULL_ALPHAA( CODES, CATEGORIES, CRITERION )
% Calculate Aickin's Alpha for nominal coding from two raters
%   [ALPHAA] = FULL_ALPHAA(CODES, CATEGORIES, CRITERION)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., rater).
%   This function can only handle two raters and nominal categories.
%
%   CATEGORIES is an optional parameter specifying the possible categories
%   as a numerical vector. If this variable is not specified, then the
%   possible categories are inferred from the CODES matrix. This can
%   underestimate reliability if all possible categories aren't used.
%
%   CRITERION is an optional parameter specifying the convergence criterion
%   for the iterative algorithm. The algorithm will stop when consecutive
%   alphas differ by less than this value. By default, it is set to 1e-7.
%
%   ALPHAA is a chance-adjusted index of agreement.
%
%   Example usage: [ALPHAA] = FULL_ALPHAA(smiledata,[0,1],0.001);
%
%   (c) Jeffrey M Girard, 2016
%   
%   References:
%
%   Aickin, M. (1990). Maximum likelihood estimation of agreement in the
%   constant predictive probability model, and its relation to Cohen's
%   kappa. Biometrics, 46, 293–302.
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
x(~isfinite(x)) = [];
if nargin < 2
    CATEGORIES = x;
    CRITERION = 1e-7;
elseif nargin < 3
    CRITERION = 1e-7;
end
if isempty(CATEGORIES)
    CATEGORIES = x;
end
CATEGORIES = sort(unique(CATEGORIES(:)));
q = length(CATEGORIES);
%% Output basic descriptives
fprintf('Number of items = %d\n',n);
fprintf('Number of raters = %d\n',r);
fprintf('Possible categories = %s\n',mat2str(CATEGORIES));
fprintf('Observed categories = %s\n',mat2str(x));
fprintf('Convergence criterion = %d\n',CRITERION);
%% Check for valid data
if n < 1
    ALPHAA = NaN;
    fprintf('ERROR: At least 1 item is required.\n')
    return;
end
if r ~= 2
    ALPHAA = NaN;
    fprintf('ERROR: Exactly 2 raters are required.\n')
    return;
end
if any(ismember(x,CATEGORIES)==0)
    fprintf('ERROR: Categories were observed in CODES that were not included in CATEGORIES.\n');
    return;
end
%% Construct contingency table
ct = nan(q,q);
for k = 1:q
    for l = 1:q
        ct(l,k) = sum(CODES(:,1)==k & CODES(:,2)==l);
    end
end
ct = ct + 1/numel(ct);
d = eye(q);
n = n + 1;
%% Calculate starting values
p_o = sum(sum(d .* ct)) / n;
n_kplus = sum(ct,2);
n_plusk = sum(ct,1);
p_kplus = n_kplus / n;
p_plusk = n_plusk / n;
p_c = 0;
for k = 1:q
    for l = 1:q
        p_c = p_c + p_kplus(k) * p_plusk(l) * d(k,l);
    end
end
alpha = (p_o - p_c) / (1 - p_c);
flag = 1;
%% Iterate until converence
while flag==1
    old_alpha = alpha;
    p_kplus = n_kplus ./ (n * (1 - alpha + alpha * d * p_plusk' / p_c));
    p_kplus(1) = 1 - sum(p_kplus(2:end));
    p_plusk = n_plusk ./ (n * (1 - alpha + alpha * p_kplus' * d / p_c));
    p_plusk(1) = 1 - sum(p_plusk(2:end));
    p_c = 0;
    for k = 1:q
        for l = 1:q
            p_c = p_c + p_kplus(k) * p_plusk(l) * d(k,l);
        end
    end
    alpha = (p_o - p_c) / (1 - p_c);
    flag = abs(alpha - old_alpha) > CRITERION;
end
fprintf('Percent observed agreement = %.3f\n',p_o);
fprintf('Percent chance agreement = %.3f\n',p_c);
fprintf('\nAickin''s alpha coefficient = %.3f\n',alpha);
ALPHAA = alpha;

end