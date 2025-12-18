#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import xgboost as xgb
from sklearn.metrics import r2_score,make_scorer,mean_squared_error
from sklearn.model_selection import GroupKFold
from scipy.stats import wilcoxon
from scipy import signal
import numpy as np
import os
from sklearn.base import clone
from matplotlib import pyplot as plt


# make a boosted tree model
model_base = xgb.XGBRegressor(
    objective='reg:squarederror', # Poisson loss
    n_estimators=200,
    max_depth=10,
    learning_rate=0.1,
    subsample=0.9,
    colsample_bytree=0.8,
    min_child_weight=10,
    reg_alpha=0.1,      # L1 regularization
    reg_lambda=5,     # L2 regularization
    gamma=1,          # Prune weak splits
    random_state=42,
    tree_method="hist"
)


fd_save_res = 'FD_SAVE_RES'
run_name = 'RUN_NAME'
n_fold_split = N_FOLD_SPLIT


## read data
fn_predictor = os.path.join(fd_save_res, f'{run_name}.predictor.csv')
X = np.loadtxt(fn_predictor, delimiter=',')
fn_output = os.path.join(fd_save_res, f'{run_name}.output.csv')
y = np.loadtxt(fn_output, delimiter=',')
print(X.shape, y.shape)


fn_group = os.path.join(fd_save_res, f'{run_name}.group.csv')
grp = pd.read_csv(fn_group)
groups = np.squeeze(np.asarray(grp))
print(groups.shape)


## split data into k-folds based on group ID, then train the XGBoost model, gather metrics
# data from the same ID don't get split into train and test, to prevent leakage
deltas = []
r2s_model = []
pct_improvements_full = []

gkf = GroupKFold(n_splits=n_fold_split)
for k, (tr, te) in enumerate(gkf.split(X, y, groups)):
    print('Running Fold '+str(k))
    X_tr = X[tr]
    X_te = X[te]
    y_tr = y[tr]
    y_te = y[te]
    
    # null = train-fold mean rate
    y_mean = np.mean(y_tr)
    y_null = np.full_like(y_te, y_mean, dtype=float)

    # ===== fit model =====
    model = clone(model_base)
    model.fit(X_tr, y_tr)
    y_hat = np.clip(model.predict(X_te), 1e-9, None)

    # save metrics
    r2s_model.append(r2_score(y_te, y_hat))
    s_model = -mean_squared_error(y_te, y_hat)
    s_null  = -mean_squared_error(y_te, y_null)
    deltas.append(s_model - s_null)
    pct_improvements_full.append(100.0 * (s_null - s_model) / s_null)
    
    # save predictions
    pred_res = np.array([y_te, y_hat])
    fn_pred_res = os.path.join(fd_save_res, f'{run_name}.pred_actual.fold{k}.csv')
    np.savetxt(fn_pred_res, pred_res)
    
    # also save the model for future use
    fn_model = os.path.join(fd_save_res, f'{run_name}.model.fold{k}.json')
    model.save_model(fn_model)

    
# save the metrics as dataframe
metrics = pd.DataFrame()
metrics.loc[:, 'ki'] = range(len(deltas))
metrics.loc[:, 'delta'] = deltas
metrics.loc[:, 'r2s_model'] = r2s_model
metrics.loc[:, 'improvement_over_null'] = pct_improvements_full
fn_metrics = os.path.join(fd_save_res, f'{run_name}.metrics.csv')
metrics.to_csv(fn_metrics)





