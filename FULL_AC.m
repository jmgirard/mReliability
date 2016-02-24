function [AC, P_O, P_C, SE, CI] = FULL_AC(CODES, CATEGORIES, SCALE, RATIO)
% Calculate the generalized form of Gwet's Agreement Coefficient (AC1/AC2)
%   [AC, P_O, P_C, SE, CI] = FULL_AC(CODES, CATEGORIES, SCALE, RATIO)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., coder).
%   This function can handle any number of coders and values.
%
%   CATEGORIES is an optional parameter specifying the possible categories
%   as a numerical vector. If this variable is not specified, then the
%   possible categories are inferred from the CODES matrix. This can
%   underestimate reliability if all possible categories aren't used.
%
%   SCALE is an optional parameter specifying the scale of measurement:
%   -Use 'nominal' for unordered categories (default)
%   -Use 'ordinal' for ordered categories of unequal size
%   -Use 'interval' for ordered categories with equal spacing
%   -Use 'ratio' for ordered categories with equal spacing and a zero point
%
%   RATIO is an optional parameter that can be used to specify the sampling
%   fraction for the current reliability experiment; it is used in the
%   calculation of SE and CI. To generalize from a sample of n items (this
%   is the number of items in CODES) to a population of N items (this is
%   the number of items in your entire dataset), set RATIO to the fraction
%   of n to N (i.e., n/N). The default of 0 can be used when N is unknown.
%
%   AC is a chance-adjusted index of agreement. When the nominal scale is
%   selected, Gwet calls it AC1; with other scales it is called AC2. AC
%   ranges from -1 to 1 where 0 means coders were no better than chance.
%
%   P_O is the percent observed agreement (from 0.000 to 1.000).
%
%   P_C is the estimated percent chance agreement (from 0.000 to 1.000).
%   
%   SE is the standard error, conditional on rater sample.
%
%   CI is a two-element vector containing the lower and upper bounds of
%   the 95% confidence interval for the AC estimate (based on the SE).
%
%   Example usage: [AC,P_O,P_C,SE,CI] = FULL_AC(smiledata,[0,1],'nominal',0);
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

%% Remove items with all missing codes
CODES(all(~isfinite(CODES),2),:) = [];
%% Calculate basic descriptives
[n,j] = size(CODES);
nprime = sum(sum(isfinite(CODES),2)>=2);
x = unique(CODES);
x(~isfinite(x)) = [];
if isempty(CATEGORIES)
    CATEGORIES = x;
end
CATEGORIES = sort(unique(CATEGORIES(:)));
q = length(CATEGORIES);
if nargin < 2
    SCALE = 'nominal';
    RATIO = 0;
elseif nargin < 3
    SCALE = 'nominal';
    RATIO = 0;
elseif nargin < 4
    RATIO = 0;
end
%% Output basic descriptives
fprintf('Number of items = %d\n',n);
fprintf('Number of coders = %d\n',j);
fprintf('Possible categories = %s\n',mat2str(CATEGORIES));
fprintf('Observed categories = %s\n',mat2str(x));
fprintf('Scale of measurement = %s\n',SCALE);
fprintf('Sampling fraction = %.3f\n',RATIO);
%% Check for valid data from two coders
if n < 1
    AC = NaN;
    fprintf('\nERROR: At least 1 item is required.\n')
    return;
end
if j < 2
    AC = NaN;
    fprintf('\nERROR: At least 2 coders are required.\n');
    return;
end
if any(ismember(x,CATEGORIES)==0)
    fprintf('ERROR: Categories were observed in CODES that were not included in CATEGORIES.\n');
    return;
end
%% Calculate weights
w = nan(length(x));
for k = 1:length(x)
    for l = 1:length(x)
        switch SCALE
            case 'nominal'
                w = eye(q);
            case 'ordinal'
                if k==l
                    w(k,l) = 1;
                else
                    M_kl = nchoosek((max(k,l) - min(k,l) + 1),2);
                    M_1q = nchoosek((max(1,q) - min(1,q) + 1),2);
                    w(k,l) = 1 - (M_kl / M_1q);
                end
            case 'interval'
                if k==l
                    w(k,l) = 1;
                else
                    dist = abs(x(k) - x(l));
                    maxdist = max(x) - min(x);
                    w(k,l) = 1 - (dist / maxdist);
                end
            case 'ratio'
                w(k,l) = 1 - (((x(k) - x(l)) / (x(k) + x(l)))^2) / (((max(x) - min(x)) / (max(x) + min(x)))^2);
                if x(k)==0 && x(l)==0, w(k,l) = 1; end
            otherwise
                error('Scale must be nominal, ordinal, interval, or ratio');
        end
    end
end
%% Calculate pi for each category
pi_k = zeros(length(x),1);
for k = 1:length(x)
    for i = 1:n
        r_ik = sum(CODES(i,:)==x(k));
        r_i = sum(isfinite(CODES(i,:)));
        pi_k(k) = pi_k(k) + r_ik / r_i;
    end
    pi_k(k) = pi_k(k) / n;
end
%% Calculate percent agreement for each item and overall
p_o_i = zeros(n,1);
for i = 1:n
    r_i = sum(isfinite(CODES(i,:)));
    if r_i >= 2
        for k = 1:length(x)
            r_ik = sum(CODES(i,:)==x(k));
            rstar_ik = 0;
            for l = 1:length(x)
                w_kl = w(k,l);
                r_il = sum(CODES(i,:)==x(l));
                rstar_ik = rstar_ik + (w_kl * r_il);
            end
            p_o_i(i) = p_o_i(i) + (r_ik * (rstar_ik - 1)) / (r_i * (r_i - 1));
        end
    end
end
P_O = sum(p_o_i) / sum(sum(isfinite(CODES),2) >= 2);
%% Calculate percent chance agreement for each item and overall
p_c_i = zeros(n,1);
for i = 1:n
    r_i = sum(isfinite(CODES(i,:)));
    for k = 1:length(x)
        r_ik = sum(CODES(i,:)==x(k));
        p_c_i(i) = p_c_i(i) + (r_ik / r_i) * (1 - pi_k(k));
    end
    T_w = sum(sum(w));
    p_c_i(i) = p_c_i(i) * (T_w / (q * (q - 1)));
end
P_C = mean(p_c_i);
%% Calculate AC point estimate
AC = (P_O - P_C) / (1 - P_C);
%% Calculate the variance of the AC point estimate
v_inner = nan(n,1);
for i = 1:n
    r_i = sum(isfinite(CODES(i,:)));
    if r_i < 2
        gamma_i = 0;
    else
        gamma_i = (n / nprime) * (p_o_i(i) - P_C) / (1 - P_C);
    end
    gammastar_i = gamma_i - 2 * (1 - AC) * ((p_c_i(i) - P_C) / (1 - P_C));
    v_inner(i) = (gammastar_i - AC) ^ 2;
end
v = ((1 - RATIO) / n) * (1 / (n - 1)) * sum(v_inner);
%% Calculate the standard error and confidence interval
SE = sqrt(v);
CI = [AC - 1.96 * SE, AC + 1.96 * SE];
%% Output reliability and variance components
fprintf('Percent observed agreement = %.3f\n',P_O);
fprintf('Percent chance agreement = %.3f\n',P_C);
fprintf('\nGwet''s AC = %.3f\n',AC);
fprintf('Standard Error (SE) = %.3f\n',SE);
fprintf('95%% Confidence Interval = %.3f to %.3f\n',CI(1),CI(2));

end