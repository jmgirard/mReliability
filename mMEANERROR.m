function [RESULTS] = mMEANERROR(DATA, YRANGE, YMEAN, YSD)
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
%   YSD is an optional parameter representing the standard deviation (SD)
%   of the measure to be used in normalizing MAE and RMSE. This can be
%   helpful if the typical SD is known from prior work. If not specified,
%   the observed SD of the second column of DATA will be used.
%
%   RESULTS is a struct containing measures derived from DATA:
%       -MAE: Mean Absolute Error
%       -NMAE_R: MAE Normalized by Range
%       -NMAE_M: MAE Normalized by Mean
%       -NMAE_S: MAE Normalized by Standard Deviation
%       -RMSE: Root Mean Square Error
%       -NRMSE_R: RMSE Normalized by Range
%       -NRMSE_M: RMSE Normalized by Mean
%       -NRMSE_S: RMSE Normalized by Standard Deviation
%
%   Example usage: mMEANERROR([preds, labels])
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

if nargin < 4
    YSD = nanstd(CRITERION);
else
    if ~isfinite(YSD)
        error('YSD must be a finite number');
    end
end

RESULTS.MAE = nanmean(abs(CRITERION - TEST));
RESULTS.NMAE_R = RESULTS.MAE / YRANGE;
RESULTS.NMAE_M = RESULTS.MAE / YMEAN;
RESULTS.NMAE_S = RESULTS.MAE / YSD;

RESULTS.RMSE = sqrt(nanmean((CRITERION - TEST) .^ 2));
RESULTS.NRMSE_R = RESULTS.RMSE / YRANGE;
RESULTS.NRMSE_M = RESULTS.RMSE / YMEAN;
RESULTS.NRMSE_S = RESULTS.RMSE / YSD;

end