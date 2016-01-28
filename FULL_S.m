function [S, P_A, P_E, SE, CI] = FULL_S(DATA, SCALE, RATIO)
% Calculate the generalized form of the S index and its properties
%   [S, P_A, P_E, SE, CI] = FULL_S(DATA, TYPE, RATIO)
%
%	DATA is a numerical matrix of codes with missing data indicated by NaN.
%   Each row is a single item and each column is a single coder. Any number
%   of coders and categories can be used. Codes should be made numerical
%   (e.g., LOW=>1, MEDIUM=>2, HIGH=>3; APPLE=>1, ORANGE=>2, BANANA=>3).
%
%   SCALE is a string corresponding to the type of weights to be used.
%   'nominal' weights are used for unordered categories
%   'ordinal' weights are used for ordered categories of unequal size
%   'interval' weights are used for ordered categories with equal spacing
%   'ratio' weights are used for ordered categories with a meaningful zero
%
%	RATIO is the sampling fraction for the current reliability experiment.
%	To generalize from a sample of n items to a population of N items, set
%   RATIO to the fraction of n to N (i.e., n/N). If the population size is
%   unknown, then use the default of 0. This parameter is optional.
%
%	S is a chance-corrected agreement index that estimates chance
%	agreement using a category-based approach.
%
%	P_A is the percent agreement observed between coders.
%
%	P_E is the percent agreement expected to be due to chance. 
%   
%	SE is the standard error of the estimate (conditional on rater sample).
%
%	CI is a two-element vector containing the lower and upper bounds of
%	the 95% confidence interval for the S estimate (based on SE).
%
%   (c) Jeffrey M Girard, 2015
%   
%   Reference: Gwet, K. L. (2014). Handbook of inter-rater reliability:
%   The definitive guide to measuring the extent of agreement among raters
%   (4th ed.). Gaithersburg, MD: Advanced Analytics.

%% Calculate variables
if nargin < 3, RATIO = 0; end
DATA(all(~isfinite(DATA),2),:) = [];
[n,~] = size(DATA);
x = unique(DATA(:));
x(~isfinite(x)) = [];
q = length(x);

%% Calculate weights
w = nan(q,q);
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
                    dist = abs(x(k) - x(l));
                    maxdist = max(x) - min(x);
                    w(k,l) = 1 - (dist / maxdist);
                end
            case 'ratio'
                w(k,l) = 1 - (((x(k) - x(l)) / (x(k) + x(l)))^2) / (((max(x) - min(x)) / (max(x) + min(x)))^2);
                if x(k)==0 && x(l)==0, w(k,l) = 1; end
            otherwise
                error('Type must be nominal, ordinal, interval, or ratio');
        end
    end
end

%% Calculate percent agreement for each item and overall
p_a_i = zeros(n,1);
for i = 1:n
    r_i = sum(isfinite(DATA(i,:)));
    if r_i >= 2
        for k = 1:q
            r_ik = sum(DATA(i,:)==x(k));
            rstar_ik = 0;
            for l = 1:q
                w_kl = w(k,l);
                r_il = sum(DATA(i,:)==x(l));
                rstar_ik = rstar_ik + (w_kl * r_il);
            end
            p_a_i(i) = p_a_i(i) + (r_ik * (rstar_ik - 1)) / (r_i * (r_i - 1));
        end
    end
end
P_A = sum(p_a_i) / sum(sum(isfinite(DATA),2)>=2);

%% Calculate percent chance agreement for each item and overall
p_e_i = zeros(n,1);
for i = 1:n
    T_w = sum(sum(w));
    p_e_i(i) = T_w / (q ^ 2);
end
P_E = mean(p_e_i);

%% Calculate K point estimate
S = (P_A - P_E) / (1 - P_E);

%% Return if variance is not requested
if nargout <=3
    SE = NaN;
    CI = [NaN,NaN];
    return;
end

%% Calculate variance of S point estimate
S_i = nan(n,1);
v_inner = 0;
for i = 1:n
    r_i = sum(isfinite(DATA(i,:)));
    if r_i >= 2
        nprime = sum(~isnan(p_a_i));
        S_i(i) = (n / nprime) * (p_a_i(i) - P_E) / (1 - P_E);
    else
        S_i(i) = 0;
    end
    v_inner = v_inner + (S_i(i) - S) ^ 2;
end
v = ((1 - RATIO) / n) * (1 / (n - 1)) * sum(v_inner);

%% Calculate the standard error and confidence interval
SE = sqrt(v);
CI = [S - 1.96 * SE, S + 1.96 * SE];

end