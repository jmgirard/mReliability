function [ ALPHAA ] = FAST_ALPHAA( CODES, CRITERIA )
% Calculate Aickin's Alpha for nominal coding from two raters
%   [ALPHAA] = FAST_ALPHAA(CODES)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., rater).
%   This function can only handle two raters and nominal categories.
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
%   Gwet, K. L. (2014). Handbook of inter-rater reliability: The definitive
%   guide to measuring the extent of agreement among raters (4th ed.).
%   Gaithersburg, MD: Advanced Analytics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate initial variables
[n,r] = size(CODES);
q = length(unique(CODES(:)));
%% Output basic descriptives
fprintf('Number of items = %d\n',n);
fprintf('Number of raters = %d\n',r);
fprintf('Number of possible categories = %d\n',2);
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
%% Construct contingency table
ctable = nan(q,q);
for k = 1:q
    for l = 1:q
        ctable(l,k) = sum(CODES(:,1)==k & CODES(:,2)==l);
    end
end
%% Initiate values
alpha = nan(1001,1);
p_k_A = nan(1001,q);
p_k_B = nan(1001,q);
p_c = nan(1001,1);
%% Calculate constants
p_o = trace(ctable) / n;
p_kplus = sum(ctable,1) ./ n;
p_plusk = (sum(ctable,2) ./ n)';
%% Calculate starting values
p_k_A(1,:) = p_kplus;
p_k_B(1,:) = p_plusk;
p_c(1) = p_k_A(1,:) * p_k_B(1,:)';
alpha(1) = (p_o - p_c(1)) / (1 - p_c(1));
%% Iterate until convergence
for t = 1:1000
    p_c(t) = p_k_A(t,:) * p_k_B(t,:)';
    for k = 1:q
        p_k_A(t+1,k) = p_kplus(k) / ((1 - alpha(t)) + (alpha(t) * p_k_B(t,k))/(p_c(t)));
        p_k_B(t+1,k) = p_plusk(k) / ((1 - alpha(t)) + (alpha(t) * p_k_A(t,k))/(p_c(t)));
    end
    alpha(t+1) = (p_o - p_c(t)) / (1 - p_c(t));
    if t==1
        continue;
    else
        if abs(alpha(t) - alpha(t+1)) < CRITERIA
            ALPHAA = alpha(t+1);
            fprintf('Converged in %d iterations\n',t);
            fprintf('Percent observed agreement = %.3f\n',p_o);
            fprintf('Percent chance agreement = %.3f\n',p_c(t));
            fprintf('\nAickin''s alpha coefficient = %.3f\n',alpha(t+1));
            return;
        end
    end
end
fprintf('Failed to converge in 1000 iterations\n');
fprintf('Percent observed agreement = %.3f\n',p_o);
fprintf('Final percent chance agreement = %.3f\n',p_c(t));
fprintf('\nFinal Aickin''s alpha coefficient = %.3f\n',alpha(t+1));

end