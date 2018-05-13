function [RESULTS] = mCONTINGENCY(DATA, TYPE, POSCLASS, NEGCLASS)
% Construct a classification contingency table and many derived measures
%
%   DATA is a numeric matrix that either contains a 2-by-2 contingency
%   table in the format of [TP, FP; FN, TN] or the raw prediction data such
%   that each row is one object, the first column contains the test
%   results, and the second column contains the trusted criterion labels.
%
%   TYPE is a required parameter indicating what the DATA variable is:
%       'table' indicates that DATA contains a 2-by-2 contingency table
%       'raw' indicates that DATA contains n-by-2 prediction data
%
%   POSCLASS is an optional parameter for the 'raw' TYPE indicating the
%   value in DATA that corresponds to the positive class (default = 1).
%
%   NEGCLASS is an optional parameter for the 'raw' TYPE indicating the
%   value in DATA that corresponds to the negative class (default = 0).
%
%   RESULTS is a struct containing measures derived from DATA:
%       -Table: a 2x2 double matrix containing [TP, FP; FN, TN]
%       -ACC: Accuracy, Agreement
%       -TPR: True Positive Rate, Sensitivity, Recall, Hit Rate
%       -TNR: True Negative Rate, Specificity
%       -PPV: Positive Predictive Value, Precision
%       -NPV: Negative Predictive Value
%       -FNR: False Negative Rate, Miss Rate, (1 - TPR)
%       -FPR: False Positive Rate, Fall-out, (1 - TNR)
%       -FDR: False Discovery Rate, (1 - PPV)
%       -FOR: False Omission Rate, (1 - NPV)
%       -PLR: Positive Likelihood Rate
%       -NLR: Negative Likelihood Rate
%       -F1S: F1 Score
%       -MCC: Matthews Correlation Coefficient, Phi Coefficient
%       -BMI: Bookmaker Informedness
%       -MRK: Markedness
%       -DOR: Diagnostic Odds Ratio
%
%   Example usage: mCONTINGENCY([50, 10; 20, 20], 'table')
%   Example usage: mCONTINGENCY([preds, labels], 'raw', 1, 0)
%
%   (c) Jeffrey M. Girard, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isnumeric(DATA)
    error('data must be numeric');
end

if nargin < 2
    error('type is required');
end

if nargin < 3
    POSCLASS = 1;
    NEGCLASS = 0;
end

if strcmpi(TYPE, 'table')
    TP = DATA(1, 1);
    FP = DATA(1, 2);
    FN = DATA(2, 1);
    TN = DATA(2, 2);
elseif strcmpi(TYPE, 'raw')
    test = DATA(:, 1);
    criterion = DATA(:, 2);
    TP = sum((test == POSCLASS) & (criterion == POSCLASS));
    FP = sum((test == POSCLASS) & (criterion == NEGCLASS));
    FN = sum((test == NEGCLASS) & (criterion == POSCLASS));
    TN = sum((test == NEGCLASS) & (criterion == NEGCLASS));
else
    error('type must be ''table'' or ''raw''');
end

RESULTS.Table = [TP, FP; FN, TN];

RESULTS.ACC = (TP + TN) / (TP + TN + FP + FN);
RESULTS.TPR = TP / (TP + FN);
RESULTS.TNR = TN / (TN + FP);
RESULTS.PPV = TP / (TP + FP);
RESULTS.NPV = TN / (TN + FN);
RESULTS.FNR = 1 - RESULTS.TPR;
RESULTS.FPR = 1 - RESULTS.TNR;
RESULTS.FDR = 1 - RESULTS.PPV;
RESULTS.FOR = 1 - RESULTS.NPV;
RESULTS.PLR = RESULTS.TPR / RESULTS.FPR;
RESULTS.NLR = RESULTS.FNR / RESULTS.TNR;

RESULTS.F1S = (2 * TP) / (2 * TP + FP + FN);
RESULTS.MCC = (TP * TN - FP * FN) / ...
    sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
RESULTS.BMI = RESULTS.TPR + RESULTS.TNR - 1;
RESULTS.MRK = RESULTS.PPV + RESULTS.NPV - 1;
RESULTS.DOR = RESULTS.PLR / RESULTS.NLR;

end