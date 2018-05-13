function [RESULTS] = mMEANERROR(DATA, YRANGE, YMEAN)
% Calculate versions of mean absolute error and root mean square error
%
%   DATA is a numeric matrix that contains the raw prediction data such
%   that each row is one object, the first column contains the test
%   results, and the second column contains the trusted criterion labels.
%
%   YRANGE is an optional parameter representing the range of the measure
%   to be used in normalizing MAE and RMSE. This can be helpful if the
%   range is restricted in principle or if the typical range is known from
%   prior work. If not specified, the observed range of the second column
%   of DATA will be used.
%
%   YMEAN is an optional parameter representing the mean of the measure to
%   be used in normalizing MAE and RMSE. This can be helpful if the typical
%   mean is known from prior work. If not specified, the observed mean of
%   the second column of DATA will be used.
%
%   RESULTS is a struct containing measures derived from DATA:
%       -MAE: Mean Absolute Error
%       -NMAE: MAE Normalized by Range
%       -CVMAE: MAE Normalized by Mean
%       -RMSE: Root Mean Square Error
%       -NRMSE: RMSE Normalized by Range
%       -CVRMSE: RMSE Normalized by Mean
%
%   Example usage: mMEANERROR([preds, labels], 100, 0)
%
%   (c) Jeffrey M. Girard, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isnumeric(DATA)
    error('DATA must be a numeric matrix')
end

if size(DATA, 2) ~= 2
    error('DATA must have 2 columns')
end

TEST = DATA(:, 1);
CRITERION = DATA(:, 2);

if nargin < 2
    YRANGE = range(CRITERION);
else
    if ~isfinite(YRANGE)
        error('YRANGE must be a finite number');
    end
end

if nargin < 3
    YMEAN = nanmean(CRITERION);
else
    if ~isfinite(YMEAN)
        error('YMEAN must be a finite number');
    end
end

RESULTS.MAE = nanmean(abs(CRITERION - TEST));
RESULTS.NMAE = RESULTS.MAE / YRANGE;
RESULTS.CVMAE = RESULTS.MAE / YMEAN;

RESULTS.RMSE = sqrt(nanmean((CRITERION - TEST) ^ 2));
RESULTS.NRMSE = RESULTS.RMSE / YRANGE;
RESULTS.CVRMSE = RESULTS.RMSE / YMEAN;

end