function [K, p_a, p_e, SE, CI, IMP, BM] = KAPPAQ(M, type, f)
% Calculate the Bennett/Brennan-Prediger Agreement Coefficient
%   [K, p_a, p_e, SE, CI, IMP, BM] = KAPPAQ(M, type, f)
%
%   (c) Jeffrey M Girard, 2015
%   
%	M is a numerical matrix of measurements. Each row is an object of
%	measurement and each column is a single source of measurement.
%	Measurements should be made numerical (e.g., A=>1, B=>2, C=>3).
%	Missing measurements should be indicated by NaN, Inf, or -Inf.
%
%   type is a string corresponding to the type of weights to be used.
%   'nominal' weights are used for unordered categories
%   'ordinal' weights are used for ordered categories of unequal size
%   'linear' weights are used for ordered categories with equal spacing
%   'quadratic' weights are used for ordered categories with equal spacing
%   'radical' weights are used for ordered categories with equal spacing
%   'ratio' weights are used for ordered categories with a meaningful zero
%
%	f is the sampling fraction for the current reliability experiment.
%	To generalize from a sample of n objects to a population of N objects,
%	set f equal to the fraction of n to N (i.e., f = n / N). If the size
%	of the population of objects is unknown, then use the default of 0.
%
%	K is an index of the agreement between two or more sources of
%	measurement that corrects for agreement due to chance.
%
%	p_a is the percent agreement observed between sources of measurement.
%
%	p_e is the percent agreement expected to be due to random guessing. 
%   
%	SE is the standard error of the K estimate (conditional on rater sample).
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
%   Reference: Bennett, E. M., Alpert, R., & Goldstein, A. C. (1954).
%   Communication through limited response questioning.
%   The Public Opinion Quarterly, 18(3), 303–308.

%% Calculate variables
if nargin < 3, f = 0; end
M(all(~isfinite(M),2),:) = [];
[n,~] = size(M);
x = unique(M(:));
x(~isfinite(x)) = [];
q = length(x);

%% Calculate weights
w = nan(q,q);
for k = 1:q
    for l = 1:q
        switch type
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
    r_i = sum(isfinite(M(i,:)));
    if r_i >= 2
        for k = 1:q
            r_ik = sum(M(i,:)==x(k));
            rstar_ik = 0;
            for l = 1:q
                w_kl = w(k,l);
                r_il = sum(M(i,:)==x(l));
                rstar_ik = rstar_ik + (w_kl * r_il);
            end
            p_a_i(i) = p_a_i(i) + (r_ik * (rstar_ik - 1)) / (r_i * (r_i - 1));
        end
    end
end
p_a = sum(p_a_i) / sum(sum(isfinite(M),2)>=2);

%% Calculate percent chance agreement for each item and overall
p_e_i = zeros(n,1);
for i = 1:n
    T_w = sum(sum(w));
    p_e_i(i) = T_w / (q ^ 2);
end
p_e = mean(p_e_i);

%% Calculate K point estimate
K = (p_a - p_e) / (1 - p_e);

%% Return if variance is not requested
if nargout <=3
    SE = NaN; CI = [NaN,NaN]; IMP = NaN; BM = '';
    return;
end

%% Calculate variance of K point estimate
K_i = nan(n,1);
v_inner = 0;
for i = 1:n
    r_i = sum(isfinite(M(i,:)));
    if r_i >= 2
        nprime = sum(~isnan(p_a_i));
        K_i(i) = (n / nprime) * (p_a_i(i) - p_e) / (1 - p_e);
    else
        K_i(i) = 0;
    end
    v_inner = v_inner + (K_i(i) - K) ^ 2;
end
v = ((1 - f) / n) * (1 / (n - 1)) * sum(v_inner);

%% Calculate the standard error and confidence interval
SE = sqrt(v);
CI = [K - 1.96 * SE, K + 1.96 * SE];

%% Calculate the interval membership probabilities
intervals = [0.8, 1.0; 0.6, 0.8; 0.4, 0.6; 0.2, 0.4; -1, 0.2];
IMP = nan(1,size(intervals,1));
for z = 1:size(intervals,1)
    a = intervals(z,1);
    b = intervals(z,2);
    IMP(z) = normcdf((K - a) / SE) - normcdf((K - b) / SE);
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