function [RESULT] = mABSERROR(DATA, NORMALIZER)
% Calculate mean absolute error with different normalization options
%
%   DATA is a numeric matrix where each row is one object, the first column
%   contains the test results (predictions), and the second column contains
%   the criterion results (trusted labels). Note that this function is most
%   appropriate when the data being compared is dimensional.
%
%   NORMALIZER is an optional parameter that specifies what the mean
%   absolute error value should be divided by (default = 1). Valid options
%   include 'range', 'mean', 'std', and any real number. If 'range',
%   'mean', or 'std' is chosen, the values are calculated from the
%   observed criterion results (i.e., the second column of DATA).
%
%   RESULT is a scalar containing the mean absolute error if NORMALIZER
%   was set to 1, or the normalized mean absolute error otherwise.
%
%   Example usage: mABSERROR([preds, labels])
%   Example usage: mABSERROR([preds, labels], 'range')
%   Example usage: mABSERROR(readydata, 50.23)
%
%   (c) Jeffrey M. Girard, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isnumeric(DATA)
    disp('DATA must be a numeric matrix');
    RESULT = NaN;
    return;
end

if size(DATA, 2) ~= 2
    disp('DATA must have 2 columns');
    RESULT = NaN;
    return;
end

if nargin < 2
    NORMALIZER = 1;
end

TEST = DATA(:, 1);
CRIT = DATA(:, 2);

if strcmpi(NORMALIZER, 'range')
    NORMALIZER = range(CRIT);
elseif strcmpi(NORMALIZER, 'mean')
    NORMALIZER = nanmean(CRIT);
elseif strcmpi(NORMALIZER, 'std')
    NORMALIZER = nanstd(CRIT);
end

MAE = nanmean(abs(CRIT - TEST));
RESULT = MAE / NORMALIZER;

end