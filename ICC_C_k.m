function [ICC, LB, UB] = ICC_C_k(DATA, ALPHA)
% Calculate the average rater consistency intraclass correlation coefficient
%   [ICC, LB, UB] = ICC_C_k(DATA)
%
%   DATA is a numerical matrix of ratings (missing values = NaN).
%   Each row is a single item and each column is a single rater.
%
%   ALPHA is the Type I error rate for the confidence interval (optional).
%
%   ICC is the reliability of the mean of the ratings from the included
%   raters. Reliability is gauged as consistency (not absolute agreement).
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
MSR = max([0, tbl{3, 4}]);
MSE = max([0, tbl{4, 4}]);
%% Calculate average rater consistency ICC
[n, k] = size(DATA);
ICC = (MSR - MSE) / MSR;
%% Calculate the confidence interval if requested
if nargout > 1
    if nargin < 2
        ALPHA = 0.05;
    end
    FL = (MSR / MSE) / finv((1 - ALPHA / 2), (n - 1), (n - 1) * (k - 1));
    FU = (MSR / MSE) / finv((1 - ALPHA / 2), (n - 1) * (k - 1), (n - 1));
    LB = 1 - 1 / FL;
    UB = 1 - 1 / FU;
end

end
