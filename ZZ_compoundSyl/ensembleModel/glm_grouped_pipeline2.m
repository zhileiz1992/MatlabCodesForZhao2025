function [mdl, metrics, coefTbl, testPredTbl] = glm_grouped_pipeline2(d_plot, pred_names, output_name, group_name, dist, link, train_test_ratio, offset_name)
%GLM_GROUPED_PIPELINE  Group-wise split GLM with held-out evaluation.
%
% Inputs
%   d_plot           : table with predictors, response, group column, and (optionally) offset
%   pred_names       : 1×P cell array of predictor names
%   output_name      : response variable name
%   dist             : distribution (e.g., 'normal','gamma','poisson', ...)
%   link             : link function (e.g., 'identity','log', ...)
%   train_test_ratio : scalar in (0,1), fraction of GROUPS assigned to training
%   group_name       : column name holding group IDs (same group must stay in one split)
%   offset_name      : (optional) name of exposure variable for log-offset; pass "" or [] if none
%
% Outputs
%   mdl         : fitted GeneralizedLinearModel on training data
%   metrics     : struct with test-set metrics:
%                   - General: MSE, R2_SSE (1 - SSE/SST), corr, corrR2, nTrain/nTest, groups
%                   - If Poisson: Dev_model, Dev_null, R2_dev_test, meanDeviance, meanLogLik, dispersion_test
%   coefTbl     : model coefficients table (Estimate, SE, tStat, pValue, etc.)
%   testPredTbl : test rows with [group, predictors..., actual, predicted]
%
% Notes
% - Rows with missing values in any used variable are removed before splitting.
% - For Poisson with offset, we use log(offset) both in fitting and prediction.
% - Test deviance explained (R2_dev_test) and mean predictive log-likelihood
%   are appropriate for Poisson GLMs and reported when dist == 'poisson' and link == 'log'.
%
% Example (Poisson counts with bin width 'exposure'):
%   rng(1);
%   [mdl,metrics,coefTbl,te] = glm_grouped_pipeline( ...
%       d_plot, {'t_onset','vae1','vae2'}, 'spk', 'poisson','log', 0.8, 'sessionID', 'exposure');

    arguments
        d_plot table
        pred_names (1,:) cell
        output_name {mustBeTextScalar}
        group_name {mustBeTextScalar}
        dist {mustBeTextScalar} = 'normal'
        link {mustBeTextScalar} = 'identity'
        train_test_ratio (1,1) double {mustBeGreaterThan(train_test_ratio,0), mustBeLessThan(train_test_ratio,1)} = 0.8
        offset_name = ""   % string/char or empty
    end

    % --- 1) Keep only needed columns; drop rows with missing
    usedVars = [pred_names, {output_name}, {group_name}];
    if ~isempty(offset_name) && any(strcmpi(string(offset_name), d_plot.Properties.VariableNames))
        usedVars = [usedVars, {char(offset_name)}];
        hasOffset = true;
    else
        hasOffset = false;
    end
    T = d_plot(:, usedVars);
    T = rmmissing(T);

    if height(T) == 0
        error('No rows remain after removing missing values in required variables.');
    end

    % --- 2) Group-wise split
    G = T.(group_name);
    uG = unique(G);
    nG = numel(uG);
    idxPerm = randperm(nG);
    nTrainG = max(1, min(nG-1, round(train_test_ratio * nG)));
    trainGroups = uG(idxPerm(1:nTrainG));
    isTrain = ismember(G, trainGroups);
    Ttr = T(isTrain, :);
    Tte = T(~isTrain, :);

    if isempty(Tte)
        warning('Test set is empty after split; metrics will be NaN.');
    end

    % --- 3) Fit GLM on training data
    formula = sprintf('%s ~ %s', output_name, strjoin(pred_names, ' + '));
    if hasOffset
        mdl = fitglm(Ttr, formula, 'Distribution', char(dist), 'Link', char(link), ...
                     'Offset', log(Ttr.(offset_name)));
    else
        mdl = fitglm(Ttr, formula, 'Distribution', char(dist), 'Link', char(link));
    end

    % --- 4) Predict on test set
    if ~isempty(Tte)
        y_true = Tte.(output_name);
        if hasOffset
            y_pred = predict(mdl, Tte, 'Offset', log(Tte.(offset_name)));
        else
            y_pred = predict(mdl, Tte);
        end
    else
        y_true = [];
        y_pred = [];
    end

    % --- 5) Metrics (general + Poisson-specific when applicable)
    metrics = struct();
    metrics.nTest          = numel(y_true);
    metrics.nTrain         = height(Ttr);
    metrics.nGroupsTrain   = nTrainG;
    metrics.nGroupsTotal   = nG;

    % General (SSE/SST-style; interpret cautiously for non-normal models)
    if ~isempty(y_true)
        residuals = y_true - y_pred;
        SSE = sum(residuals.^2);
        SST = sum( (y_true - mean(y_true)).^2 );
        metrics.MSE    = mean(residuals.^2);
        metrics.R2_SSE = 1 - SSE / max(eps, SST);
        if numel(y_true) > 1
            metrics.corr   = corr(y_true, y_pred, 'Rows','complete');
            metrics.corrR2 = metrics.corr^2;
        else
            metrics.corr   = NaN; metrics.corrR2 = NaN;
        end
    else
        metrics.MSE = NaN; metrics.R2_SSE = NaN; metrics.corr = NaN; metrics.corrR2 = NaN;
    end

    % Poisson/log diagnostics (held-out)
    if ~isempty(y_true) && strcmpi(dist,'poisson') && strcmpi(link,'log')
        y = y_true;
        mu = max(y_pred, eps);

        % Null model mean on TEST set
        if hasOffset
            expo = Tte.(offset_name);
            expo = max(expo, eps);
            rate0 = sum(y) / sum(expo);         % MLE for intercept-only with offset
            mu0   = rate0 * expo;                % mu0 = exp(b0)*exposure
        else
            mu0 = repmat(mean(y), size(y));     % intercept-only without offset
        end
        mu0 = max(mu0, eps);

        % Poisson deviance on test
        % (Define 0*log(0/.) as 0 via max(y,eps))
        Dev_model = 2 * sum( y .* log( max(y,eps)./mu ) - (y - mu) );
        Dev_null  = 2 * sum( y .* log( max(y,eps)./mu0) - (y - mu0) );

        % Predictive log-likelihood (per observation)
        meanLogLik   = mean( y .* log(mu) - mu - gammaln(y + 1) );
        meanDeviance = Dev_model / numel(y);

        % A dispersion heuristic on test (uses model DF from training)
        dfe = max(1, numel(y) - mdl.NumEstimatedCoefficients);
        dispersion_test = Dev_model / dfe;

        metrics.Dev_model      = Dev_model;
        metrics.Dev_null       = Dev_null;
        metrics.R2_dev_test    = 1 - Dev_model / max(eps, Dev_null);
        metrics.meanDeviance   = meanDeviance;
        metrics.meanLogLik     = meanLogLik;
        metrics.dispersion_test = dispersion_test;
    else
        metrics.Dev_model = NaN; metrics.Dev_null = NaN;
        metrics.R2_dev_test = NaN; metrics.meanDeviance = NaN;
        metrics.meanLogLik = NaN; metrics.dispersion_test = NaN;
    end

    % --- 6) Coefficients and significance
    coefTbl = mdl.Coefficients;

    % --- 7) Per-row outputs for the test set
    if ~isempty(Tte)
        testPredTbl = Tte(:, [{group_name}, pred_names, {output_name}]);
        testPredTbl.Properties.VariableNames{end} = 'actual';
        testPredTbl.predicted = y_pred;
    else
        testPredTbl = Tte; % empty table
    end
end


