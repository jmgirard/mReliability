function [RESULT] = mDES(GROUP1, GROUP2, TYPE)
% Calculate the standardized mean difference effect size
%
%   GROUP1 and GROUP2 are numeric vectors from separate samples, groups, 
%   time points, or contexts that you want to compare across. If TYPE is
%   'diff' then GROUP1 and GROUP2 must have the same length.
%
%   TYPE is an optional parameter that specifies what standardizer should
%   be used in the denominator of the standardized mean difference.
%       'pool': The pooled SD using Cohen's formula (default)
%       's1': The SD observed in GROUP1
%       's2': The SD observed in GROUP2
%       'diff': The SD of the observed differences between groups
%
%   RESULT is a scalar containing the standardized mean difference (d).
%
%   Example usage: mDES(females, males)
%   Example usage: mDES(time1, time2, 'diff')
%   Example usage: mDES(patients, controls, 's2')
%
%   (c) Jeffrey M. Girard, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isvector(GROUP1) || ~isvector(GROUP2) || ~isnumeric(GROUP1) || ~isnumeric(GROUP2)
    disp('GROUP1 and GROUP2 must be numeric vectors.');
    RESULT = NaN;
    return;
end

if nargin < 3
    TYPE = 'pool';
end

n1 = nansum(GROUP1);
n2 = nansum(GROUP2);
m1 = nanmean(GROUP1);
m2 = nanmean(GROUP2);
s1 = nanstd(GROUP1);
s2 = nanstd(GROUP2);
contrast = m1 - m2;

if strcmpi(TYPE, 'pool')
    v1 = nanvar(GROUP1);
    v2 = nanvar(GROUP2);
    vpool = (v1 * (n1 - 1) + v2 * (n2 - 1)) / (n1 + n2 - 2);
    spool = sqrt(vpool);
    RESULT = contrast / spool;
elseif strcmpi(TYPE, 's1')
    RESULT = contrast / s1;
elseif strcmpi(TYPE, 's2')
    RESULT = contrast / s2;
elseif strcmpi(TYPE, 'diff')
    if length(GROUP1) == length(GROUP2)
        mdiff = nanmean(GROUP1 - GROUP2);
        sdiff = nanstd(GROUP1 - GROUP2);
        RESULT = mdiff / sdiff;
    else
        RESULT = NaN;
    end
else
    disp('TYPE must be ''pool'', ''s1'', ''s2'', or ''diff''.');
    RESULT = NaN;
end

end