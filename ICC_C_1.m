function [ICC,LB,UB] = ICC_C_1(DATA,ALPHA)
% Calculate the single rater consistency intraclass correlation coefficient
%   [ICC,LB,UB] = ICC_C_1(DATA)
%
%   DATA is a numerical matrix of ratings with no missing values.
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
%% Calculate the mean sums of squares and ICC(C,1)
MSR = SSR / (n - 1);
MSE = SSE / ((n - 1)*(k - 1));
ICC = (MSR - MSE) / (MSR + MSE*(k-1));
%% Calculate the confidence interval if requested
if nargout > 1
    if nargin < 2
        ALPHA = 0.05;
    end
    FL = (MSR/MSE) / finv((1-ALPHA/2),(n - 1),(n - 1)*(k - 1));
    FU = (MSR/MSE) / finv((1-ALPHA/2),(n - 1)*(k - 1),(n - 1));
    LB = (FL - 1) / (FL + k - 1);
    UB = (FU - 1) / (FU + k - 1);
end

end