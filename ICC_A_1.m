function [ICC,LB,UB] = ICC_A_1(DATA,ALPHA)
% Calculate the single rater agreement intraclass correlation coefficient
%   [ICC,LB,UB] = ICC_A_1(DATA)
%
%   DATA is a numerical matrix of ratings with no missing values.
%   Each row is a single item and each column is a single rater.
%   This coefficient is appropriate for interval or ratio measurements.
%
%   ALPHA is the Type I error rate for the confidence interval (optional).
%
%   ICC is the reliability of the ratings taken from any single included
%   rater. Reliability is gauged as agreement on an absolute scale.
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
SSC = 0;
SSE = 0;
for i = 1:n
    for j = 1:k
        SSR = SSR + (y_i(i) - y)^2;
        SSC = SSC + (y_j(j) - y)^2;
        SSE = SSE + (DATA(i,j) - y_j(j) - y_i(i) + y)^2;
    end
end
%% Calculate the mean sums of squares and ICC(A,1)
MSR = SSR / (n - 1);
MSE = SSE / ((n - 1)*(k - 1));
MSC = SSC / (k - 1);
ICC = (MSR - MSE) / (MSR + MSE*(k - 1) + (k/n)*(MSC - MSE));
%% Calculate the confidence interval if requested
if nargout > 1
    if nargin < 2
        ALPHA = 0.05;
    end
    a   = (k*ICC) / (n*(1 - ICC));
    b   = 1 + (k*ICC*(n - 1))/(n*(1 - ICC));
    v   = ((a*MSC + b*MSE)^2) / (((a*MSC)^2)/(k - 1) + ((b*MSE)^2)/((n - 1)*(k - 1)));
    FL = finv(1 - ALPHA/2,n - 1,v);
    FU = finv(1 - ALPHA/2,v,n - 1);
    LB  = (n*(MSR - FL*MSE)) / (FL*(k*MSC + MSE*(k*n - k - n)) + n*MSR);
    UB  = (n*(FU*MSR - MSE)) / (k*MSC + MSE*(k*n - k - n) + n*FU*MSR);
end

end