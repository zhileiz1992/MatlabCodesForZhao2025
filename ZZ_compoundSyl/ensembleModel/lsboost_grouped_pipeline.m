function [mdl, metrics, importTbl, testPredTbl] = lsboost_grouped_pipeline( ...
        d_plot, pred_names, output_name, group_name, train_test_ratio)
%LSBOOST_GROUPED_PIPELINE  Group-wise split + LSBoost regression (fitrensemble).
%
% Inputs
%   d_plot           : table with predictors, response, and group column
%   pred_names       : 1×P cell array of predictor variable names
%   output_name      : response variable name (numeric, continuous)
%   train_test_ratio : scalar in (0,1), fraction of GROUPS assigned to training
%   group_name       : column name holding group IDs (same group must stay in one split)
%
% Outputs
%   mdl         : RegressionEnsemble (LSBoost) trained on the training split
%   metrics     : struct with test-set MSE, R2_SSE, corr, corrR2, sizes
%   importTbl   : table of predictor importances (descending)
%   testPredTbl : test rows with [group, predictors..., actual, predicted]
%
% Notes
% - Set RNG before calling for reproducible splits, e.g., rng(1).
% - Categorical predictors are supported if columns are categorical in d_plot.

    arguments
        d_plot table
        pred_names (1,:) cell
        output_name {mustBeTextScalar}
        group_name {mustBeTextScalar}
        train_test_ratio (1,1) double {mustBeGreaterThan(train_test_ratio,0), mustBeLessThan(train_test_ratio,1)} = 0.8
    end

    % --- 1) Keep only needed columns; drop rows with missing
    usedVars = [pred_names, {output_name}, {group_name}];
    T = d_plot(:, usedVars);
    T = rmmissing(T);
    if height(T) == 0
        error('No rows remain after removing missing values in required variables.');
    end

    % --- 2) Group-wise split to avoid leakage
    G  = T.(group_name);
    uG = unique(G);
    nG = numel(uG);
    idxPerm  = randperm(nG);
    nTrainG  = max(1, min(nG-1, round(train_test_ratio * nG)));
    trainGs  = uG(idxPerm(1:nTrainG));
    isTrain  = ismember(G, trainGs);
    Ttr = T(isTrain, :);
    Tte = T(~isTrain, :);

    if isempty(Tte)
        warning('Test set is empty after split; metrics will be NaN.');
    end

    % --- 3) Fit LSBoost (squared-error boosting of regression trees)
    nCycles   = 300;                                % boosting iterations
    learnRate = 0.1;                                % shrinkage
%     maxSplits = max(2, round(0.5 * numel(pred_names) + 10));  % depth proxy
    maxSplits = 10;
    minLeaf   = max(1, floor(height(Ttr)/200));

    tTree = templateTree('MaxNumSplits', maxSplits, ...
                         'MinLeafSize',  minLeaf);

    mdl = fitrensemble( ...
            Ttr, output_name, ...
            'Method','LSBoost', ...
            'PredictorNames', pred_names, ...
            'Learners', tTree, ...
            'NumLearningCycles', nCycles, ...
            'LearnRate', learnRate);

    % --- 4) Predict on test set
    if ~isempty(Tte)
        y_true = Tte.(output_name);
        y_pred = predict(mdl, Tte);
    else
        y_true = []; y_pred = [];
    end

    % --- 5) Test-set metrics
    metrics = struct('nTest', numel(y_true), ...
                     'nTrain', height(Ttr), ...
                     'nGroupsTrain', nTrainG, ...
                     'nGroupsTotal', nG, ...
                     'MSE', NaN, 'R2_SSE', NaN, ...
                     'corr', NaN, 'corrR2', NaN);
    if ~isempty(y_true)
        resid = y_true - y_pred;
        SSE = sum(resid.^2);
        SST = sum((y_true - mean(y_true)).^2);
        metrics.MSE    = mean(resid.^2);
        metrics.R2_SSE = 1 - SSE / max(eps, SST);
        if numel(y_true) > 1
            metrics.corr   = corr(y_true, y_pred, 'Rows','complete');
            metrics.corrR2 = metrics.corr.^2;
        end
    end

    % --- 6) Predictor importances
    try
        imp = predictorImportance(mdl);   % 1×P
        importTbl = table(pred_names(:), imp(:), ...
                          'VariableNames', {'Predictor','Importance'});
        importTbl = sortrows(importTbl, 'Importance', 'descend');
    catch
        importTbl = table(pred_names(:), NaN(numel(pred_names),1), ...
                          'VariableNames', {'Predictor','Importance'});
    end

    % --- 7) Per-row outputs for the test set
    if ~isempty(Tte)
        testPredTbl = Tte(:, [{group_name}, pred_names, {output_name}]);
        testPredTbl.Properties.VariableNames{end} = 'actual';
        testPredTbl.predicted = y_pred;
    else
        testPredTbl = Tte;  % empty
    end
end
