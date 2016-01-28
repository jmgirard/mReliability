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
if size(CODES,2) ~= 2
    error('FAST_S can only handle two columns.');
end
v = unique(CODES);
if length(v) ~= 2
    error('FAST_S can only handle two categories.');
end
%% Calculate contingency table cells
a = sum(CODES(:,1)==v(1) & CODES(:,2)==v(1));
b = sum(CODES(:,1)==v(2) & CODES(:,2)==v(1));
c = sum(CODES(:,1)==v(1) & CODES(:,2)==v(2));
d = sum(CODES(:,1)==v(2) & CODES(:,2)==v(2));
%% Calculate observed and chance agreement
p_o = (a + d) / (a + b + c + d);
p_c = 1 / 2;
%% Calculate reliability index
S = (p_o - p_c) / (1 - p_c);

end

