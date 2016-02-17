function [S] = FAST_S(CODES, Q)
% Quickly calculate the S index using nominal data from exactly two coders
%   [S] = FAST_S(CODES, Q)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., coder).
%   This function requires nominal coding from exactly two coders (the
%   FULL_S function can handle any number of coders and different scales).
%   Items with any missing codes (e.g., NaN) will be dropped from analysis.
%
%   Q is an optional parameter specifying the total number of categories.
%   If this variable is not specified, then the number of categories is
%   inferred from the CODES matrix. This inference can underestimate S if
%   all possible categories aren't included in CODES.
%
%   S is a chance-adjusted index of agreement. It assumes that each
%   category has an equal chance of being selected at random. It ranges
%   from -1.0* to 1.0 where 0.0 means coders were no better than chance.
%   *The actual lower bound is determined by the number of categories.
%   
%   Example usage: S = FAST_S(smiledata,2);
%   
%   (c) Jeffrey M Girard, 2016
%
%   References:
%   
%   Bennett, E. M., Alpert, R., & Goldstein, A. C. (1954).
%   Communication through limited response questioning.
%   The Public Opinion Quarterly, 18(3), 303–308.
%
%   Zhao, X., Liu, J. S., & Deng, K. (2012).
%   Assumptions behind inter-coder reliability indices.
%   In C. T. Salmon (Ed.), Communication Yearbook (pp. 418–480). Routledge.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remove items with any missing codes
CODES(any(~isfinite(CODES),2),:) = [];
%% Calculate basic descriptives
[n,j] = size(CODES);
x = unique(CODES);
if nargin < 2
    q = length(x);
else
    q = Q;
end
%% Output basic descriptives
fprintf('Number of items = %d\n',n);
fprintf('Number of coders = %d\n',j);
fprintf('Number of possible categories = %d\n',q);
fprintf('Observed categories = %s\n',mat2str(x));
%% Check for valid data from two coders
if n < 1
    S = NaN;
    fprintf('ERROR: At least 1 item is required.\n')
    return;
end
if j ~= 2
    S = NaN;
    fprintf('ERROR: Use FULL_S() for more than 2 coders.\n');
    return;
end
%% Calculate components and reliability
agr = sum(CODES(:,1)==CODES(:,2));
dis = sum(CODES(:,1)~=CODES(:,2));
p_o = agr / (agr + dis);
p_c = 1 / q;
S = (p_o - p_c) / (1 - p_c);
%% Output components and reliability
fprintf('Number of agreements = %d\n',agr);
fprintf('Number of disagreements = %d\n',dis);
fprintf('Percent observed agreement = %.3f\n',p_o);
fprintf('Percent chance agreement = %.3f\n',p_c);
fprintf('\nS index = %.3f\n',S);

end