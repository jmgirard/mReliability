function [icc, sigma_r_sq, sigma_e_sq] = ICC1B1(DATA)
% Calculate the Single Rater Intraclass Correlation under Model 1B
%   [icc, sigma_r_sq, sigma_e_sq] = ICC1B1(DATA)
%   
%   DATA is a numerical matrix where each row represents a subject and each
%   column represents an arbitrary rater position (the specific raters that
%   occupy each position may differ between subjects). Missing ratings
%   should be represented as NaN.
%   
%   icc is the single rater ICC calculated under model 1B. It is used as
%   a measure of intra-rater reliability when each rater may rate a
%   (partially or completely) different group of subjects.
%
%   sigma_r_sq is the rater variance.
%
%   sigma_e_sq is the error variance. This estimate combines variance from
%   both the subject factor and unspecified experimental error factor.
%
%   (c) Jeffrey M Girard, 2015

n = size(DATA,1);
r = size(DATA,2);
m_ij = isfinite(DATA);
m_i = sum(m_ij,2);
m_j = sum(m_ij,1);
M = sum(m_i);
y_j = nansum(DATA,1);
y_j_sq = y_j .^ 2;
T_2r = sum(y_j_sq ./ m_j);
T_y = sum(y_j);
T_2y = sum(nansum(DATA .^ 2, 2));
k_1 = 0;
for i = 1:n
    for j = 1:r
        k_1 = k_1 + (m_ij(i,j) .^ 2) / (m_i(i));
    end
end
sigma_e_sq = (T_2y - T_2r) / (M - r);
sigma_r_sq = (T_2r - (T_y ^ 2 / M) - (r - 1) * sigma_e_sq) / (M - k_1);
icc = sigma_r_sq / (sigma_r_sq + sigma_e_sq);

end