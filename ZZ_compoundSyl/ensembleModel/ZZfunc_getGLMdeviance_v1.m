function [R2_dev_test] = ZZfunc_getGLMdeviance_v1(mdl, X_te, y_te)
% evaluate trained GLM model (mdl) on test datasets
% for Poisson distribution with log link only

% make prediction
y_pred = predict(mdl, X_te);

% calculate deviance
y = table2array(y_te);
mu = max(y_pred, eps);
mu0 = repmat(mean(y), size(y));
mu0 = max(mu0, eps);

% Poisson deviance on test
% (Define 0*log(0/.) as 0 via max(y,eps))
Dev_model = 2 * sum( y .* log( max(y,eps)./mu ) - (y - mu) );
Dev_null  = 2 * sum( y .* log( max(y,eps)./mu0) - (y - mu0) );
R2_dev_test = 1 - Dev_model / max(eps, Dev_null);

end

