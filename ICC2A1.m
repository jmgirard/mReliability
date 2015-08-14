function [icc, icc_a, sigma_s_sq, sigma_r_sq, sigma_sr_sq, sigma_e_sq] = ICC2A1(DATA)
% Calculate Single Rater Intraclass Correlations under Model 2
%   [icc, icc_a, sigma_s_sq, sigma_r_sq, sigma_sr_sq, sigma_e_sq] = ICC2A1(DATA)
%   
%   DATA is a numerical matrix where each row represents a subject and each
%   column represents an arbitrary rater position (the specific raters that
%   occupy each position may differ between subjects). The first column
%   should contain integers corresponding to subject IDs. This allows some
%   subjects to be rated multiple times by one or more raters. Because of 
%   this possibility, the subject-rater interaction is included. Missing
%   ratings are permitted and should be represented in DATA as NaN.
%   
%   icc is the single rater inter-rater reliability ICC calculated under
%   model 2 with the subject-rater interaction included.
%
%   icc_a is the intra-rater reliability ICC calculated under model 2.
%
%   sigma_s_sq is the subject variance.
%
%   sigma_r_sq is the rater variance.
%
%   sigma_sr_sq is the subject-rater interaction variance.
%
%   sigma_e_sq is the error variance. This estimate combines variance from
%   both the subject factor and unspecified experimental error factor.
%
%   (c) Jeffrey M Girard, 2015

x_ijk = DATA(:,1);
x_val = unique(x_ijk);
y_ijk = DATA(:,2:end);
n = length(x_val);
r = size(y_ijk,2);

y_ij = nan(n,r);
m_ij = nan(n,r);
for i = 1:n
    x = x_val(i);
    for j = 1:r
        index = x_ijk==x;
        y_ij(i,j) = nansum(y_ijk(index,j));
        m_ij(i,j) = sum(isfinite(y_ijk(index,j)));
    end
end
m_i = sum(m_ij,2);
m_j = sum(m_ij,1);
M = sum(m_i);
m_i_sq = m_i .^ 2;
m_j_sq = m_j .^ 2;

k_1 = sum(m_i_sq);
k_1p = k_1 / M;
k_2 = sum(m_j_sq);
k_2p = k_2 / M;
k_3 = 0;
for i = 1:n
    k_3 = k_3 + sum((m_ij(i,:) .^ 2) ./ m_i(i));
end
k_4 = 0;
for j = 1:r
    k_4 = k_4 + sum((m_ij(:,j) .^ 2) ./ m_j(j));
end
k_5 = sumsqr(m_ij);
k_5p = k_5 / M;

y_i = sum(y_ij,2);
y_i_sq = y_i .^ 2;
y_j = sum(y_ij,1);
y_j_sq = y_j .^ 2;
T_2r = sum(y_j_sq ./ m_j);
T_2sr = sum(nansum((y_ij .^ 2) ./ m_ij));
T_2y = nansum(y_ijk(:) .^ 2);
T_2s = sum(y_i_sq ./ m_i);
lambda_0 = sum(m_ij(:) > 0);
T_y = sum(y_j);
T_yp_sq = (T_y .^ 2) / M;

sigma_e_sq = (T_2y - T_2sr) / (M - lambda_0);
delta_s = (T_2sr - T_2s - sigma_e_sq * (lambda_0 - n)) / (M - k_3);
delta_r = (T_2sr - T_2r - sigma_e_sq * (lambda_0 - r)) / (M - k_4);
sigma_sr_sq = (delta_r * (M - k_1p) + delta_s * (k_3 - k_2p) - (T_2s - T_yp_sq - sigma_e_sq * (n - 1))) / (M - k_1p - k_2p + k_5p);
sigma_r_sq = delta_s - sigma_sr_sq;
sigma_s_sq = delta_r - sigma_sr_sq;

icc = sigma_s_sq / (sigma_s_sq + sigma_r_sq + max([0,sigma_sr_sq]) + sigma_e_sq);
icc_a = (sigma_r_sq + sigma_s_sq + max([0,sigma_sr_sq])) / (sigma_s_sq + sigma_r_sq + max([0,sigma_sr_sq]) + sigma_e_sq);

end