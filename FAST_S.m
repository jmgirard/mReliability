function [S] = FAST_S(CODES,Q)
% Quickly calculate the S index using data from exactly two coders
%   [S] = FAST_S(CODES,Q)
%
%   CODES should be a numerical matrix where each row corresponds to a
%   single item of measurement (e.g., participant or question) and each
%   column corresponds to a single source of measurement (i.e., coder).
%   This function requires exactly two coders (the function FULL_S can
%   handle any number of coders). Codes can take on any number of values.
%
%   Q is an optional parameter that can be used to specify the number of
%   possible values. If this variable is not specified, then the number
%   of possible values is inferred from the CODES matrix. This inference 
%   can underestimate S if all possible values aren't included in CODES.
%
%   S is a chance-corrected index of agreement. It assumes that each
%   category has an equal chance of being selected at random. It ranges
%   from -1.0 to 1.0 where 0.0 means coders were no better than chance.
%   Note that the lower bound is determined by the number of values.
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

%% Get basic descriptives
n = size(CODES,1);
j = size(CODES,2);
v = unique(CODES);
if nargin < 2
    q = length(v);
else
    q = Q;
end
%% Output basic descriptives
fprintf('Number of items = %d\n',n);
fprintf('Number of coders = %d\n',j);
fprintf('Number of values = %d\n',q);
fprintf('Observed values = %s\n',mat2str(v));
%% Check for binary coding from two observers
if n < 1
    S = NaN;
    fprintf('S = NaN; At least 1 item is required.\n')
    return;
end
if j ~= 2
    S = NaN;
    fprintf('S = NaN; Use FULL_S() for more than 2 coders.\n');
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
fprintf('S index = %.3f\n',S);
end