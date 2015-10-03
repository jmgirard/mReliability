function [ICC,LB,UB] = ICC_C_k(DATA,ALPHA)
% Calculate the average rater consistency intraclass correlation coefficient
%   [ICC,LB,UB] = ICC_C_k(DATA)
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
[rowindex,~] = find(~isfinite(DATA));
DATA(rowindex,:) = [];
%% Calculate descriptive statistics
[n,k] = size(DATA);
y = mean(DATA(:));
y_j = mean(DATA,1);
y_i = mean(DATA,2);
%% Calculate row, column, and error sums of squares
SSR = 0;
SSE = 0;
for i = 1:n
    for j = 1:k
        SSR = SSR + (y_i(i) - y)^2;
        SSE = SSE + (DATA(i,j) - y_j(j) - y_i(i) + y)^2;
    end
end
%% Calculate the mean sums of squares and ICC(C,k)
MSR = SSR / (n - 1);
MSE = SSE / ((n - 1)*(k - 1));
ICC = (MSR - MSE) / MSR;
%% Calculate the confidence interval if requested
if nargout > 1
    if nargin < 2
        ALPHA = 0.05;
    end
    FL = (MSR/MSE) / finv((1-ALPHA/2),(n - 1),(n - 1)*(k - 1));
    FU = (MSR/MSE) / finv((1-ALPHA/2),(n - 1)*(k - 1),(n - 1));
    LB  = 1 - 1/FL;
    UB  = 1 - 1/FU;
end

end