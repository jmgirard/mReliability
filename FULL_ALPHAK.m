function [ALPHA, P_O, P_C, SE, CI] = FULL_ALPHAK(CODES, CATEGORIES, SCALE, RATIO)
% Calculate Krippendorff's alpha coefficient using generalized formulas
%   [ALPHA, P_O, P_C, SE, CI] = FULL_ALPHAK(CODES, CATEGORIES, SCALE, RATIO)
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
%   ALPHA is a chance-adjusted index of agreement.
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
%   Example usage: [ALPHA,P_O,P_C,SE,CI] = FULL_ALPHA(smiledata,[0,1],'nominal',0);
%
%   (c) Jeffrey M Girard, 2016
%   
%   References:
%   
%   Krippendorff, K. (1970). Estimating the reliability, systematic error
%   and random error of interval data. Educational and Psychological
%   Measurement, 30(1), 61–70.
%   
%   Krippendorff, K. (1980). Content analysis: An introduction to its
%   methodology. Newbury Park, CA: Sage Publications.
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
if exist(CATEGORIES,'var')==0
    CATEGORIES = x;
elseif isempty(CATEGORIES)
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
%% Check for valid data from multiple coders
if n < 1
    ALPHA = NaN;
    fprintf('\nERROR: At least 1 item is required.\n')
    return;
end
if j < 2
    ALPHA = NaN;
    fprintf('\nERROR: At least 2 coders are required.\n');
    return;
end
if any(ismember(x,CATEGORIES)==0)
    fprintf('ERROR: Categories were observed in CODES that were not included in CATEGORIES.\n');
    return;
end
%% Calculate weights
w = nan(q);
for k = 1:q
    for l = 1:q
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
                    dist = abs(CATEGORIES(k) - CATEGORIES(l));
                    maxdist = max(CATEGORIES) - min(CATEGORIES);
                    w(k,l) = 1 - (dist / maxdist);
                end
            case 'ratio'
                w(k,l) = 1 - (((CATEGORIES(k) - CATEGORIES(l)) / (CATEGORIES(k) + CATEGORIES(l)))^2) / (((max(CATEGORIES) - min(CATEGORIES)) / (max(CATEGORIES) + min(CATEGORIES)))^2);
                if CATEGORIES(k)==0 && CATEGORIES(l)==0, w(k,l) = 1; end
            otherwise
                error('Scale must be nominal, ordinal, interval, or ratio');
        end
    end
end
%% Calculate rbar
r_i = nan(1,n);
for i = 1:n
    r_i(i) = sum(isfinite(CODES(i,:)));
end
rbar = nanmean(r_i);
%% Calculate pi for each category
pi_k = zeros(q,1);
for k = 1:q
    for i = 1:n
        r_ik = sum(CODES(i,:)==CATEGORIES(k));
        r_i = sum(isfinite(CODES(i,:)));
        if r_i < 2
            continue;
        else
            pi_k(k) = pi_k(k) + r_ik / rbar;
        end
    end
    pi_k(k) = pi_k(k) / nprime;
end
%% Calculate percent agreement for each item and overall
p_o_i = zeros(n,1);
for i = 1:n
    r_i = sum(isfinite(CODES(i,:)));
    if r_i >= 2
        for k = 1:q
            r_ik = sum(CODES(i,:)==CATEGORIES(k));
            rstar_ik = 0;
            for l = 1:q
                w_kl = w(k,l);
                r_il = sum(CODES(i,:)==CATEGORIES(l));
                rstar_ik = rstar_ik + (w_kl * r_il);
            end
            p_o_i(i) = p_o_i(i) + (r_ik * (rstar_ik - 1)) / (rbar * (r_i - 1));
        end
    end
end
p_o = sum(p_o_i) / nprime;
%% Adjust percent observed agreement sample size for each item and overall
epsilon_n = 1 / (nprime * rbar);
p_oEn_i = zeros(n,1);
for i = 1:n
    r_i = sum(isfinite(CODES(i,:)));
    if r_i >= 2
        p_oEn_i(i) = (1 - epsilon_n) * (p_o_i(i) - p_o * (r_i - rbar) / rbar) + epsilon_n;
    end
end
P_O = mean(p_oEn_i);
%% Calculate percent chance agreement for each item and overall
p_c_i = zeros(n,1);
for i = 1:n
    r_i = sum(isfinite(CODES(i,:)));
    for k = 1:q
        r_ik = sum(CODES(i,:)==CATEGORIES(k));
        pibar_kdot = 0;
        pibar_dotk = 0;
        for l = 1:q
            pibar_kdot = pibar_kdot + w(k,l) * pi_k(l);
            pibar_dotk = pibar_dotk + w(k,l) * pi_k(k);
        end
        pibar_k = (pibar_kdot + pibar_dotk) / 2;
        p_c_i(i) = p_c_i(i) + pibar_k * r_ik / rbar;
    end
    p_c_i(i) = p_c_i(i) - (r_i - rbar) / rbar;
end
P_C = mean(p_c_i);
%% Calculate reliability point estimate
ALPHA = (P_O - P_C) / (1 - P_C);
%% Calculate the variance of the reliability point estimate
v_inner = nan(n,1);
for i = 1:n
    r_i = sum(isfinite(CODES(i,:)));
    if r_i < 2
        alpha_i = 0;
    else
        alpha_i = (p_oEn_i(i) - P_C) / (1 - P_C);
    end
    alphastar_i = alpha_i - (1 - ALPHA) * ((p_c_i(i) - P_C) / (1 - P_C));
    v_inner(i) = (alphastar_i - ALPHA) ^ 2;
end
v = ((1 - RATIO) / nprime) * (1 / (nprime - 1)) * sum(v_inner);
%% Calculate the standard error and confidence interval
SE = sqrt(v);
CI = [ALPHA - 1.96 * SE, ALPHA + 1.96 * SE];
%% Output reliability and variance components
fprintf('Percent observed agreement = %.3f\n',P_O);
fprintf('Percent chance agreement = %.3f\n',P_C);
fprintf('\nKrippendorff''s alpha = %.3f\n',ALPHA);
fprintf('Standard Error (SE) = %.3f\n',SE);
fprintf('95%% Confidence Interval = %.3f to %.3f\n',CI(1),CI(2));

end