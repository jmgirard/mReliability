function [icc, sigma_s_sq, sigma_e_sq] = ICC1A1(DATA)
% Calculate the Single Rater Intraclass Correlation under Model 1A
%   [icc, sigma_s_sq, sigma_e_sq] = ICC1A1(DATA)
%   
%   DATA is a numerical matrix where each row represents a subject and each
%   column represents an arbitrary rater position (the specific raters that
%   occupy each position may differ between subjects). Missing ratings
%   should be represented as NaN.
%   
%   icc is the single rater ICC calculated under model 1A. It is used as
%   a measure of inter-rater reliability when each subject is rated by
%   multiple raters but the group of raters differs between subjects.
%
%   sigma_s_sq is the subject variance.
%
%   sigma_e_sq is the error variance. This estimate combines variance from
%   both the rater factor and unspecified experimental error factor.
%   
%   (c) Jeffrey M Girard, 2015

n = size(DATA,1);
r = size(DATA,2);
m_ij = isfinite(DATA);
m_i = sum(m_ij,2);
m_j = sum(m_ij,1);
M = sum(m_i);
y_i = nansum(DATA,2);
y_i_sq = y_i .^ 2;
T_2s = sum(y_i_sq ./ m_i);
T_y = sum(y_i);
T_2y = sum(nansum(DATA .^ 2, 2));
k_0 = 0;
for i = 1:n
    for j = 1:r
        k_0 = k_0 + (m_ij(i,j) .^ 2) / (m_j(j));
    end
end
sigma_e_sq = (T_2y - T_2s) / (M - n);
sigma_s_sq = (T_2s - (T_y ^ 2 / M) - (n - 1) * sigma_e_sq) / (M - k_0);
icc = sigma_s_sq / (sigma_s_sq + sigma_e_sq);

end