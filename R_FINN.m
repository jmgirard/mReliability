function [R_1, R_k] = R_FINN(DATA)
% Calculate Finn's reliability coefficient
%   [R_1, R_k] = R_FINN(DATA)
%   
%   DATA is a numerical matrix of ratings (missing values = NaN).
%   Each row is a single item and each column is a single rater.
%   
%   R_1 is the reliability of any single rater of those included
%
%   R_k is the reliability of the average of all included raters
%   
%   (c) Jeffrey M Girard, 2015

%% Remove any missing values
[rowindex,~] = find(~isfinite(DATA));
DATA(rowindex,:) = [];
%% Calculate residual mean square from one-way ANOVA
[n,k] = size(DATA);
[~,tbl,~] = anova1(DATA',[],'off');
MSE = tbl{3,4};
%% Calculate single and average rater reliability
R_1 = 1.0 - MSE / ((k^2 - 1) / 12);
R_k = (n * R_1) / (1 + R_1 * (n - 1));

end