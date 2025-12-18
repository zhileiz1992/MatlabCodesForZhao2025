function [mdl, metrics, coefTbl, testPredTbl] = glm_grouped_pipeline( ...
        d_plot, pred_names, output_name, group_name, dist, link, train_test_ratio)
%GLM_GROUPED_PIPELINE  Train/test GLM with group-wise split (no leakage).
%
% Inputs
%   d_plot           : table with predictors, response, and group column
%   pred_names       : 1×P cell array of predictor variable names (strings/chars)
%   output_name      : response variable name (string/char)
%   dist             : distribution for fitglm (e.g., 'normal','gamma','poisson', ...)
%   link             : link function for fitglm (e.g., 'identity','log', ...)
%   train_test_ratio : scalar in (0,1), fraction of GROUPS assigned to training
%   group_name       : column name holding group IDs (same group must stay in one split)
%
% Outputs
%   mdl         : fitted GeneralizedLinearModel on training data
%   metrics     : struct with test-set MSE, R2 (variance explained), and corrR2
%   coefTbl     : model coefficients table (Estimate, SE, tStat, pValue, etc.)
%   testPredTbl : table for test set with actual vs predicted and group IDs
%
% Notes
% - Rows with missing values in any used variable are removed before splitting.
% - R2 is computed on held-out data as 1 - SSE/SST; corrR2 is squared Pearson
%   correlation between y_true and y_pred (useful for non-identity links).
%
% Example
%   [mdl,metrics,coefTbl,testPredTbl] = glm_grouped_pipeline( ...
%       d_plot, {'t_onset','vae1','vae2'}, 'ifr', 'normal','identity', 0.8, 'sessionID');

    arguments
        d_plot table
        pred_names (1,:) cell
        output_name {mustBeTextScalar}
        group_name {mustBeTextScalar}
        dist {mustBeTextScalar} = 'normal'
        link {mustBeTextScalar} = 'identity'
        train_test_ratio (1,1) double {mustBeGreaterThan(train_test_ratio,0), mustBeLessThan(train_test_ratio,1)} = 0.8
    end

    % --- 1) Clean rows with missing values in used variables ---
    usedVars = [pred_names, {output_name}, {group_name}];
    T = d_plot(:, usedVars);
    T = rmmissing(T);  % remove any row with missing among used variables

    if height(T) == 0
        error('No rows remain after removing missing values in required variables.');
    end

    % --- 2) Group-wise train/test split (no leakage) ---
    G = T.(group_name);
    uG = unique(G);
    nG = numel(uG);

    % Reproducibility left to caller (set rng before calling if desired)
    idxPerm = randperm(nG);
    nTrainG = max(1, min(nG-1, round(train_test_ratio * nG))); % ensure non-empty test if possible
    trainGroups = uG(idxPerm(1:nTrainG));

    isTrain = ismember(G, trainGroups);
    Ttr = T(isTrain, :);
    Tte = T(~isTrain, :);

    if isempty(Tte)
        warning('Test set ended up empty (all groups in train). Adjust train_test_ratio.');
    end

    % --- 3) Fit GLM on training data ---
    formula = sprintf('%s ~ %s', output_name, strjoin(pred_names, ' + '));
    mdl = fitglm(Ttr, formula, 'Distribution', char(dist), 'Link', char(link));

    % --- 4) Predict on test set ---
    if ~isempty(Tte)
        y_true = Tte.(output_name);
        y_pred = predict(mdl, Tte);
    else
        y_true = [];
        y_pred = [];
    end

    % --- 5) Metrics on test set ---
    metrics = struct('MSE', NaN, 'R2', NaN, 'corrR2', NaN, ...
                     'nTest', numel(y_true), 'nTrain', height(Ttr), ...
                     'nGroupsTrain', nTrainG, 'nGroupsTotal', nG);
    if ~isempty(y_true)
        residuals = y_true - y_pred;
        SSE = sum(residuals.^2);
        SST = sum( (y_true - mean(y_true)).^2 );
        metrics.MSE = mean(residuals.^2);
        metrics.R2  = 1 - SSE / max(eps, SST);   % variance explained on held-out
        if numel(y_true) > 1
            c = corr(y_true, y_pred, 'Rows', 'complete');
            metrics.corrR2 = c^2;
        end
    end

    % --- 6) Coefficients and significance ---
    coefTbl = mdl.Coefficients;   % includes Estimate, SE, tStat, pValue, etc.

    % --- 7) Assemble per-row outputs for test set ---
    if ~isempty(Tte)
        testPredTbl = Tte(:, [{group_name}, pred_names, {output_name}]);
        testPredTbl.predicted = y_pred;
        testPredTbl.Properties.VariableNames{end-1} = 'actual'; % rename response column to 'actual'
        % Move columns order to: group, predictors..., actual, predicted
        % (already in that order)
    else
        testPredTbl = Tte; % empty table
    end
end

