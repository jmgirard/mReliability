function [ICC, LB, UB] = ICC_C_1(DATA, ALPHA)
% Calculate the single rater consistency intraclass correlation coefficient
%   [ICC, LB, UB] = ICC_C_1(DATA)
%
%   DATA is a numerical matrix of ratings (missing values = NaN).
%   Each row is a single item and each column is a single rater.
%
%   ALPHA is the Type I error rate for the confidence interval (optional).
%
%   ICC is the reliability of the ratings taken from any single included
%   rater. Reliability is gauged as consistency (not absolute agreement).
%
%   LB and UB are the confidence interval's lower and upper bounds.
%
%   (c) Jeffrey M Girard, 2015
%
%   Reference: McGraw, K. O., & Wong, S. P. (1996).
%   Forming inferences about some intraclass correlation coefficients. 
%   Psychological Methods, 1(1), 30–46.

%% Remove any missing values
[rowindex, ~] = find(~isfinite(DATA));
DATA(rowindex, :) = [];
%% Calculate mean squares from two-way ANOVA
[~, tbl, ~] = anova2(DATA, 1, 'off');
MSR = tbl{3, 4};
MSE = tbl{4, 4};
%% Calculate single rater consistency ICC
[n, k] = size(DATA);
ICC = (MSR - MSE) / (MSR + MSE * (k - 1));
%% Calculate the confidence interval if requested
if nargout > 1
    if nargin < 2
        ALPHA = 0.05;
    end
    FL = (MSR / MSE) / finv((1 - ALPHA / 2), (n - 1), (n - 1) * (k - 1));
    FU = (MSR / MSE) / finv((1 - ALPHA / 2), (n - 1) * (k - 1), (n - 1));
    LB = (FL - 1) / (FL + k - 1);
    UB = (FU - 1) / (FU + k - 1);
end

end
