function [ ALPHAA ] = FAST_ALPHAA( CODES )
% Calculate Aickin's Alpha for binary coding from two raters
%   [ALPHAA] = FAST_ALPHAA(CODES)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., rater).
%   This function can only handle two raters and binary categories.
%
%   ALPHAA is a chance-adjusted index of agreement.
%
%   Example usage: [ALPHAA] = FAST_ALPHAA(smiledata);
%
%   (c) Jeffrey M Girard, 2016
%   
%   References:
%
%   Aickin, M. (1990). Maximum likelihood estimation of agreement in the
%   constant predictive probability model, and its relation to Cohen's
%   kappa. Biometrics, 46, 293–302.
%
%   Guggenmoos-Holzmann, I. (1996). The meaning of kappa: Probabilistic
%   concepts of reliability and validity revisited. Journal of Clinical
%   Epidemiology, 49(7), 775–782.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remove items with any missing codes
CODES(any(~isfinite(CODES),2),:) = [];
%% Calculate basic descriptives
[n,r] = size(CODES);
x = unique(CODES);
x(~isfinite(x)) = [];
q = length(x);
%% Output basic descriptives
fprintf('Number of items = %d\n',n);
fprintf('Number of raters = %d\n',r);
fprintf('Observed categories = %s\n',mat2str(x));
%% Check for valid data
if n < 1
    ALPHAA = NaN;
    fprintf('ERROR: At least 1 item is required.\n')
    return;
end
if r ~= 2
    ALPHAA = NaN;
    fprintf('ERROR: Exactly 2 raters are required.\n')
    return;
end
if q ~= 2
    ALPHAA = NaN;
    fprintf('ERROR: Exactly 2 categories are required.\n')
    return;
end
%% Calculate alpha coefficient
n_11 = sum(CODES(:,1)==x(1) & CODES(:,2)==x(1));
n_22 = sum(CODES(:,1)==x(2) & CODES(:,2)==x(2));
n_12 = sum(CODES(:,1)==x(1) & CODES(:,2)==x(2));
n_21 = sum(CODES(:,1)==x(2) & CODES(:,2)==x(1));
p_o = (n_11 + n_22) / n;
ALPHAA = p_o * (1 - 1 / sqrt((n_11 * n_22) / (n_12 * n_21)));
fprintf('Percent observed agreement = %.3f\n',p_o);
fprintf('Percent chance agreement = %.3f\n',p_c);
fprintf('\nAickin''s alpha coefficient = %.3f\n',alpha);

end