function [KAPPA, P_O, P_C, SE, CI] = FULL_KAPPA(CODES, CATEGORIES, SCALE, RATIO)
% Calculate Cohen's kappa coefficient using generalized formulas
%   [KAPPA, P_O, P_C, SE, CI] = FULL_KAPPA(CODES, CATEGORIES, SCALE, RATIO)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., rater).
%   This function can handle any number of raters and values.
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
%   KAPPA is a chance-adjusted index of agreement.
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
%   Example usage: [KAPPA,P_O,P_C,SE,CI] = FULL_KAPPA(smiledata,[0,1],'nominal',0);
%
%   (c) Jeffrey M Girard, 2016
%   
%   References:
%   
%   Cohen, J. (1960). A coefficient of agreement for nominal scales.
%   Educational and Psychological Measurement, 20(1), 37–46.
%
%   Cohen, J. (1968). Weighted kappa: Nominal scale agreement with
%   provision for scaled disagreement or partial credit. Psychological
%   Bulletin, 70(4), 213–220.
%
%   Conger, A. J. (1980). Integration and generalization of kappas for
%   multiple raters. Psychological Bulletin, 88(2), 322–328.
%   
%   Gwet, K. L. (2014). Handbook of inter-rater reliability: The definitive
%   guide to measuring the extent of agreement among raters (4th ed.).
%   Gaithersburg, MD: Advanced Analytics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remove items with all missing codes
CODES(all(~isfinite(CODES),2),:) = [];
%% Calculate basic descriptives
[n,r] = size(CODES);
nprime = sum(sum(isfinite(CODES),2)>=2);
x = unique(CODES);
x(~isfinite(x)) = [];
if nargin < 2
    CATEGORIES = x;
    SCALE = 'nominal';
    RATIO = 0;
elseif nargin < 3
    SCALE = 'nominal';
    RATIO = 0;
elseif nargin < 4
    RATIO = 0;
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
fprintf('Scale of measurement = %s\n',SCALE);
fprintf('Sampling fraction = %.3f\n',RATIO);
%% Check for valid data from multiple raters
if n < 1
    KAPPA = NaN;
    fprintf('\nERROR: At least 1 item is required.\n')
    return;
end
if r < 2
    KAPPA = NaN;
    fprintf('\nERROR: At least 2 raters are required.\n');
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
            p_o_i(i) = p_o_i(i) + (r_ik * (rstar_ik - 1)) / (r_i * (r_i - 1));
        end
    end
end
P_O = sum(p_o_i) / nprime;
%% Calculate percent chance agreement for each item and overall
p_c_i = zeros(n,1);
for i = 1:n
    p_c_i_inner = nan(r,q);
    for g = 1:r
        n_g = sum(isfinite(CODES(:,g)));
        for k = 1:q
            p_gk = sum(CODES(:,g)==x(k)) / n_g;
            lambda_gk_inner = zeros(q,1);
            for l = 1:q
                p_gl = sum(CODES(:,g)==x(l)) / n_g;
                delta_gl = CODES(i,g) == x(l);
                epsilon_g = isfinite(CODES(i,g));
                lambda_gk_inner(l) = w(k,l) * (delta_gl - p_gl * (epsilon_g - n_g / n));
            end
            lambda_gk = (n / n_g) * sum(lambda_gk_inner);
            pbar_plusk = zeros(r,1);
            for gg = 1:r
                pbar_plusk = pbar_plusk + sum(CODES(:,g)==x(k)) / n_gg;
            end
            pbar_plusk = pbar_plusk / r;
            p_c_i_inner(g,k) = lambda_gk * (r * pbar_plusk - p_gk);
        end
    end
    p_c_i(i) = 1 / (r * (r -1)) * sum(sum(p_c_i_inner));
end
P_C = mean(p_c_i);
%% Calculate reliability point estimate
KAPPA = (P_O - P_C) / (1 - P_C);
%% Calculate the variance of the reliability point estimate
v_inner = nan(n,1);
for i = 1:n
    r_i = sum(isfinite(CODES(i,:)));
    if r_i < 2
        kappa_i = 0;
    else
        kappa_i = (n / nprime) * (p_o_i(i) - P_C) / (1 - P_C);
    end
    kappastar_i = kappa_i - 2 * (1 - KAPPA) * ((p_c_i(i) - P_C) / (1 - P_C));
    v_inner(i) = (kappastar_i - KAPPA) ^ 2;
end
v = ((1 - RATIO) / n) * (1 / (n - 1)) * sum(v_inner);
%% Calculate the standard error and confidence interval
SE = sqrt(v);
CI = [KAPPA - 1.96 * SE, KAPPA + 1.96 * SE];
%% Output reliability and variance components
fprintf('Percent observed agreement = %.3f\n',P_O);
fprintf('Percent chance agreement = %.3f\n',P_C);
fprintf('\nCohen''s kappa coefficient = %.3f\n',KAPPA);
fprintf('Standard Error (SE) = %.3f\n',SE);
fprintf('95%% Confidence Interval = %.3f to %.3f\n',CI(1),CI(2));

end