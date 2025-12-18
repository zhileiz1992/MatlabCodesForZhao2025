#!/usr/bin/env python
# coding: utf-8

# evaluate performance of trained model on test dataset

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


model_path = 'MODEL_PATH'
fn_test_predictor = 'FN_TEST_PREDICTOR'
fn_test_output = 'FN_TEST_OUTPUT'
fn_save = 'FN_SAVE'

# model = clone(model_base)
model = xgb.XGBRegressor()
model.load_model(model_path)


## read data
X_te = np.loadtxt(fn_test_predictor, delimiter=',')
y_te = np.loadtxt(fn_test_output, delimiter=',')
print(X_te.shape, y_te.shape)

# null = train-fold mean rate
y_mean = np.mean(y_te)
y_null = np.full_like(y_te, y_mean, dtype=float)

# make predictions
y_hat = np.clip(model.predict(X_te), 1e-9, None)

# save metrics
r2s_model = r2_score(y_te, y_hat)
s_model = -mean_squared_error(y_te, y_hat)
s_null  = -mean_squared_error(y_te, y_null)
deltas = s_model - s_null
pct_improvements_full = 100.0 * (s_null - s_model) / s_null

# save metrics
metric = [r2s_model, s_model, s_null, deltas, pct_improvements_full]
fd_save = os.path.dirname(fn_save)
if not os.path.exists(fd_save):
    os.makedirs(fd_save)
np.savetxt(fn_save, metric)





