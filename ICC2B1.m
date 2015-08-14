function [icc, sigma_s_sq, sigma_r_sq, sigma_e_sq] = ICC2B1(DATA)
% Calculate Single Rater Intraclass Correlations under Model 2
%   [icc, sigma_s_sq, sigma_r_sq, sigma_e_sq] = ICC2B1(DATA)
%   
%   DATA is a numerical matrix where each row represents a subject and each
%   column represents an arbitrary rater position (the specific raters that
%   occupy each position may differ between subjects). Unlike Model 2A each
%   subject is only rated once by each rater and so the subject-rater
%   interaction is removed from this model. Missing ratings are permitted
%   and should be represented in DATA as NaN.
%   
%   icc is the single rater inter-rater reliability ICC calculated under
%   model 2 with the subject-rater interaction removed.
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

y_ij = DATA;
n = size(y_ij,1);
r = size(y_ij,2);

m_ij = ~isnan(y_ij);
m_i = sum(m_ij,2);
m_j = sum(m_ij,1);
M = sum(m_i);

k_1 = sumsqr(m_i);
k_1p = k_1 / M;
k_2 = sumsqr(m_j);
k_2p = k_2 / M;
k_3 = 0;
for i = 1:n
    k_3 = k_3 + sum((m_ij(i,:) .^ 2) ./ m_i(i));
end
k_4 = 0;
for j = 1:r
    k_4 = k_4 + sum((m_ij(:,j) .^ 2) ./ m_j(j));
end

y_i = nansum(y_ij,2);
y_i_sq = y_i .^ 2;
y_j = nansum(y_ij,1);
y_j_sq = y_j .^ 2;

T_2r = sum(y_j_sq ./ m_j);
T_2y = nansum(y_ij(:) .^ 2);
T_2s = nansum(y_i_sq ./ m_i);
T_2u = M * nanmean(y_ij(:)) .^ 2;

lambda_1 = (M - k_1p) / (M - k_4);
lambda_2 = (M - k_2p) / (M - k_3);

sigma_e_sq = (lambda_2 * (T_2y - T_2s) + lambda_1 * (T_2y - T_2r) - (T_2y - T_2u)) / (lambda_2 * (M - n) + lambda_1 * (M - r) - (M - 1));
sigma_r_sq = (T_2y - T_2s - sigma_e_sq * (M - n)) / (M - k_3);
sigma_s_sq = (T_2y - T_2r - sigma_e_sq * (M - r)) / (M - k_4);

icc = sigma_s_sq / (sigma_s_sq + sigma_r_sq + sigma_e_sq);

end