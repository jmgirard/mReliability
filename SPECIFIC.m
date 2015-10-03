function [P_S, VALUES, SE, CI] = SPECIFIC(DATA, ALPHA, NBOOT, PARALLEL)
% Calculate Specific Agreement for Categorical Annotations
%   [P_S, VALUES, SE, CI] = SPECIFIC(DATA, ALPHA, NBOOT, PARALLEL)
%
%   DATA is a numerical matrix of codes with missing data indicated by NaN.
%   Each row is a single item and each column is a single coder.
%   Codes should be numerical (e.g., NEGATIVE => -1, NEUTRAL => 0, POSITIVE => 1).
%
%   ALPHA is the Type I error rate for the confidence interval (optional).
%
%   NBOOT is the number of bootstrap samples to use during SE estimation.
%
%   PARALLEL can be set to 1 to run the bootstrap using the MATLAB
%   parallel computing toolbox (this is faster when NBOOT is large).
%
%   P_S is a row vector containing the specific agreement estimates
%   for each category, sorted in numerically-ascending order.
%
%   VALUES is a row vector containing the category labels in order.
%   
%   SE is a row vector containing the bootstrapped standard error of P_S.
%
%   CI is a matrix containing the confidence interval of P_S given ALPHA.
%   The first column contains lower bounds, the second column upper bounds.
%
%   (c) Jeffrey M Girard, 2015
%   
%   Reference: http://www.john-uebersax.com/stat/raw.htm

%% Calculate variables
[nITEMS,~] = size(DATA);
VALUES = unique(DATA(:));
nVALUES = length(VALUES);
%% Count category occurrences
n_jk = zeros(nVALUES,nITEMS);
for j = 1:nVALUES
    for k = 1:nITEMS
        n_jk(j,k) = sum(DATA(k,:)==VALUES(j));
    end
end
n_k = sum(n_jk,1);
%% Calculate specific agreement
S = zeros(nVALUES,1);
SP = zeros(nVALUES,1);
P_S = zeros(nVALUES,1);
for j = 1:nVALUES
    for k = 1:nITEMS
        S(j) = S(j) + n_jk(j,k)*(n_jk(j,k) - 1);
        SP(j) = SP(j) + n_jk(j,k)*(n_k(k) - 1);
    end
    P_S(j) = S(j) / SP(j);
end
%% Calculate bootstrapped confidence interval if requested
if nargout > 2
    if nargin < 2
        ALPHA = 0.05;
    end
    if nargin < 3
        NBOOT = 100000;
        PARALLEL = 0;
    end
    if isnan(ALPHA), SE = NaN; CI = NaN; return; end
    if PARALLEL==1
        p = parpool();
        opts = statset('UseParallel',true);
        bootstat = bootstrp(NBOOT,@SPECIFIC,DATA,NaN,'Options',opts);
        delete(p);
    else
        bootstat = bootstrp(NBOOT,@SPECIFIC,DATA,NaN);
    end
    SE = std(bootstat)';
    CI = nan(length(SE),2);
    for i = 1:length(SE)
        CI(i,1) = P_S(i) - SE(i) * norminv(1-ALPHA/2);
        CI(i,2) = P_S(i) + SE(i) * norminv(1-ALPHA/2);
    end
end

end