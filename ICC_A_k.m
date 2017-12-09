function [ICC,LB,UB] = ICC_A_k(DATA,ALPHA)
% Calculate the average rater agreement intraclass correlation coefficient
%   [ICC,LB,UB] = ICC_A_k(DATA)
%
%   DATA is a numerical matrix of ratings (missing values = NaN).
%   Each row is a single item and each column is a single rater.
%
%   ALPHA is the Type I error rate for the confidence interval (optional).
%
%   ICC is the reliability of the mean of the ratings from the included
%   raters. Reliability is gauged as agreement on an absolute scale.
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
MSC = tbl{2, 4};
MSR = tbl{3, 4};
MSE = tbl{4, 4};
%% Calculate average rater agreement ICC
[n, k] = size(DATA);
ICC = (MSR - MSE) / (MSR + (MSC - MSE) / n);
%% Calculate the confidence interval if requested
if nargout > 1
    if nargin < 2
        ALPHA = 0.05;
    end
    c  = ICC / (n * (1 - ICC));
    d  = 1 + (ICC * (n - 1)) / (n * (1 - ICC));
    v  = ((c * MSC + d * MSE) ^ 2) / (((c * MSC) ^ 2) / (k - 1) + ((d * MSE) ^ 2) / ((n - 1) * (k - 1)));
    FL = finv((1 - ALPHA / 2), (n - 1), v);
    FU = finv((1 - ALPHA / 2), v, (n - 1));
    LB = (n * (MSR - FL * MSE)) / (FL * (MSC - MSE) + n * MSR);
    UB = (n * (FU * MSR - MSE)) / (MSC - MSE + n * FU * MSR);
end

end
