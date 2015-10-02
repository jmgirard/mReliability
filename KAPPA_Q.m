function [KAPPA, P_A, P_E, SE, CI, IMP, BM] = KAPPA_Q(DATA, TYPE, RATIO)
% Calculate Kappa Q or the Free Marginal Agreement Coefficient
%   [KAPPA, P_A, P_E, SE, CI, IMP, BM] = KAPPA_Q(DATA, TYPE, RATIO)
%
%	DATA is a numerical matrix of codes with missing data indicated by NaN.
%   Each row is a single item and each column is a single coder.
%	Codes should be made numerical (e.g., APPLE=>1, ORANGE=>2, BANANA=>3).
%
%   TYPE is a string corresponding to the type of weights to be used.
%   'nominal' weights are used for unordered categories
%   'ordinal' weights are used for ordered categories of unequal size
%   'linear' weights are used for ordered categories with equal spacing
%   'quadratic' weights are used for ordered categories with equal spacing
%   'radical' weights are used for ordered categories with equal spacing
%   'ratio' weights are used for ordered categories with a meaningful zero
%
%	RATIO is the sampling fraction for the current reliability experiment.
%	To generalize from a sample of n items to a population of N items, set
%   RATIO to the fraction of n to N (i.e., n/N). If the population size is
%   unknown, then use the default of 0. This parameter is optional.
%
%	KAPPA is a chance-corrected agreement index that estimates chance
%	agreement using a category-based approach.
%
%	P_A is the percent agreement observed between coders.
%
%	P_E is the percent agreement expected to be due to chance. 
%   
%	SE is the standard error of the estimate (conditional on rater sample).
%
%	CI is a two-element vector containing the lower and upper bounds of
%	the 95% confidence interval for the K estimate.
%
%	IMP is the interval membership probability for Altman's benchmarks:
%	'Very Good':	0.8 to 1.0
%	'Good':         0.6 to 0.8
%	'Moderate':     0.4 to 0.6
%	'Fair':         0.2 to 0.4
%	'Poor':         less than 0.2
%
%	BM is the final benchmark identified by the 95% Rule, which stipulates
%	using the highest level which has a cumulative probability over 0.95
%
%   (c) Jeffrey M Girard, 2015
%   
%   Reference: Brennan, R. L., & Prediger, D. J. (1981).
%   Coefficient Kappa: Some uses, misuses, and alternatives.
%   Educational and Psychological Measurement, 41(3), 687–699.

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
        switch TYPE
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
            case 'linear'
                if k==l
                    w(k,l) = 1;
                else
                    dist = abs(x(k) - x(l));
                    maxdist = max(x) - min(x);
                    w(k,l) = 1 - (dist / maxdist);
                end
            case 'quadratic'
                if k==l
                    w(k,l) = 1;
                else
                    w(k,l) = 1 - (x(k) - x(l))^2 / (max(x) - min(x))^2;
                end
            case 'radical'
                if k==l
                    w(k,l) = 1;
                else
                    w(k,l) = 1 - sqrt(abs(x(k) - x(l))) / sqrt(max(x) - min(x));
                end
            case 'ratio'
                w(k,l) = 1 - (((x(k) - x(l)) / (x(k) + x(l)))^2) / (((max(x) - min(x)) / (max(x) + min(x)))^2);
                if x(k)==0 && x(l)==0, w(k,l) = 1; end
            otherwise
                error('Type must be nominal, ordinal, linear, quadratic, radical, or ratio');
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
KAPPA = (P_A - P_E) / (1 - P_E);

%% Return if variance is not requested
if nargout <=3
    SE = NaN; CI = [NaN,NaN]; IMP = NaN; BM = '';
    return;
end

%% Calculate variance of K point estimate
K_i = nan(n,1);
v_inner = 0;
for i = 1:n
    r_i = sum(isfinite(DATA(i,:)));
    if r_i >= 2
        nprime = sum(~isnan(p_a_i));
        K_i(i) = (n / nprime) * (p_a_i(i) - P_E) / (1 - P_E);
    else
        K_i(i) = 0;
    end
    v_inner = v_inner + (K_i(i) - KAPPA) ^ 2;
end
v = ((1 - RATIO) / n) * (1 / (n - 1)) * sum(v_inner);

%% Calculate the standard error and confidence interval
SE = sqrt(v);
CI = [KAPPA - 1.96 * SE, KAPPA + 1.96 * SE];

%% Calculate the interval membership probabilities
intervals = [0.8, 1.0; 0.6, 0.8; 0.4, 0.6; 0.2, 0.4; -1, 0.2];
IMP = nan(1,size(intervals,1));
for z = 1:size(intervals,1)
    a = intervals(z,1);
    b = intervals(z,2);
    IMP(z) = normcdf((KAPPA - a) / SE) - normcdf((KAPPA - b) / SE);
end

%% Identify the final benchmark level using the 95% rule
labels = {'Very Good','Good','Moderate','Fair','Poor'};
index = find(cumsum(IMP)>.95,1,'first');
if isempty(index)
    BM = 'Unknown';
else
    BM = labels{index};
end

end