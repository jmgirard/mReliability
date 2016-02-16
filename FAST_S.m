function [S] = FAST_S(CODES)
% Quickly calculate the S index for dichotomous coding from two coders
%   [S] = FAST_S(CODES)
%
%   CODES is a numerical matrix of codes where each row is a single item
%   and each column is a single coder. There must be two coders. Codes can
%   take on any two values, but there must be two different values present.
%
%   S is a chance-corrected index of agreement. It assumes that each
%   category has an equal chance of being selected at random. It ranges
%   from -1.0 to 1.0 where 0.0 means coders were no better than chance.
%   
%   (c) Jeffrey M Girard, 2016
%
%   Reference: Bennett, E. M., Alpert, R., & Goldstein, A. C. (1954).
%   Communication through limited response questioning.
%   The Public Opinion Quarterly, 18(3), 303–308.

%% Check for binary coding from two coders
n = size(CODES,1);
fprintf('Number of items = %d\n',n);
if n < 1
    S = NaN;
    fprintf('S = NaN; At least 1 item is required.\n')
    return;
end
j = size(CODES,2);
fprintf('Number of coders = %d\n',j);
if j ~= 2
    S = NaN;
    fprintf('S = NaN; Use FULL_S() for more than 2 coders.\n');
    return;
end
v = unique(CODES);
q = length(v);
fprintf('Number of values = %d\n',q);
fprintf('Values = %s\n',mat2str(v));
if q ~= 2
    S = NaN;
    fprintf('S = NaN; Use FULL_S() for more than 2 values.\n');
    return;
end
%% Calculate agreements and disagreements
A = sum(CODES(:,1)==CODES(:,2));
D = sum(CODES(:,1)==CODES(:,2));
fprintf('Number of agreements = %d\nNumber of disagreements = %d\n',A,D);
%% Calculate observed and chance agreement
p_o = A / (A + D);
p_c = 1 / q;
fprintf('Percent observed agreement = %.3f\nPercent chance agreement = %.3f\n',p_o,p_c);
%% Calculate reliability index
S = (p_o - p_c) / (1 - p_c);
fprintf('S index = %.3f\n',S);
end