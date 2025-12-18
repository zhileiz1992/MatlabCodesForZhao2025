# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 16:57:52 2025

@author:Glab

"""

# === System and Utility ===
import os
import glob
import warnings
import pickle
import random
from random import randint
from pathlib import Path

# === numerical packages ===
import numpy as np
import scipy.io as sio
from scipy import signal
from scipy.signal import savgol_filter, medfilt, resample_poly,firwin,filtfilt,hilbert
from scipy.integrate import cumtrapz
from scipy.interpolate import splev, splrep, interp1d
from sklearn.preprocessing import normalize
import scipy.interpolate as intp
from scipy.special import gammaln
from scipy.stats import combine_pvalues


# === PyTorch ===
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.functional import softplus, sigmoid
from torch.optim import Adam
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.utils.data import Dataset, DataLoader
from torch.distributions import MultivariateNormal, LowRankMultivariateNormal
from torchvision import transforms
from torch.utils.data import Dataset, DataLoader
from torch.utils.data import Subset, random_split

# === Plotting ===
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec

# === Loss & Activation Functions ===
SELU = nn.SELU()
mae = nn.L1Loss()

# === Data Transforms ===
p = transforms.Compose([transforms.Resize((128, 8))])

### ================= VAE ============================ # CHANGE IF thres is changed
class VAE(nn.Module):
    def __init__(self, latent_dim: int = 8, input_chan: int = 8, save_dir='', lr: float = 1e-3, model_precision: float = 10.0, device_name="auto"):
        super(VAE, self).__init__()

        self.save_dir = save_dir
        if self.save_dir != '' and not os.path.exists(self.save_dir):
            os.makedirs(self.save_dir)

        assert device_name != "cuda" or torch.cuda.is_available()
        if device_name == "auto":
            device_name = "cuda" if torch.cuda.is_available() else "cpu"
        self.device = torch.device(device_name)
        self.input_chan = input_chan
        self.latent_dim = latent_dim
        self.model_precision = model_precision
        self.lr = lr
        self._build_network()
        self.optimizer = Adam(self.parameters(), lr=self.lr)
        self.epoch = 0
        self.loss = {'train': {}, 'test': {}}
        self.to(self.device)

    def _build_network(self):

        # encoder layer
        self.conv1 = nn.Conv2d(in_channels=1, out_channels=8, kernel_size=3, stride=2, padding=2, dilation=2)
        self.conv2 = nn.Conv2d(in_channels=8, out_channels=16, kernel_size=3, stride=2, padding=2, dilation=2)
        self.conv3 = nn.Conv2d(in_channels=16, out_channels=32, kernel_size=3, stride=2, padding=2, dilation=2)

        self.bn1 = nn.BatchNorm2d(1)
        self.bn2 = nn.BatchNorm2d(8)
        self.bn3 = nn.BatchNorm2d(16)

        # MODIFIED: new flattened size is 1024 (instead of 2048) for input shape (128,16)
        self.fc1 = nn.Linear(1024, 256)  # CHANGED
        self.fc21 = nn.Linear(256, 128)
        self.fc22 = nn.Linear(256, 128)
        self.fc23 = nn.Linear(256, 128)

        self.fc31 = nn.Linear(128, self.latent_dim)
        self.fc32 = nn.Linear(128, self.latent_dim)
        self.fc33 = nn.Linear(128, self.latent_dim)

        # decoder layer
        self.fc4 = nn.Linear(self.latent_dim, 128)
        self.fc5 = nn.Linear(128, 256)
        self.fc6 = nn.Linear(256, 512)

        # MODIFIED: output to 1024 instead of 2048
        self.fc7 = nn.Linear(512, 1024)  # CHANGED

        self.convt0 = nn.ConvTranspose2d(in_channels=32, out_channels=16, kernel_size=3, stride=2, padding=1, output_padding=(1, 1))
        self.bn4 = nn.BatchNorm2d(32)
        self.convt1 = nn.ConvTranspose2d(in_channels=16, out_channels=8, kernel_size=3, stride=2, padding=1, output_padding=(1, 1))
        self.bn5 = nn.BatchNorm2d(16)
        self.convt2 = nn.ConvTranspose2d(in_channels=8, out_channels=1, kernel_size=3, stride=2, padding=1, output_padding=(1, 1))
        self.bn6 = nn.BatchNorm2d(8)

    def encoder(self, x):
        batch_size = x.size(0)

        x = x.unsqueeze(1)
        x = SELU(self.conv1(self.bn1(x)))
        x = SELU(self.conv2(self.bn2(x)))
        x = SELU(self.conv3(self.bn3(x)))

        x = x.view(batch_size, -1)
        x = SELU(self.fc1(x))

        mu = SELU(self.fc21(x))
        mu = self.fc31(mu)

        d = SELU(self.fc23(x))
        d = softplus(self.fc33(d)) + 1e-3

        return mu, d

    def decoder(self, z):

        batch_size = z.size(0)
        z = SELU(self.fc4(z))
        z = SELU(self.fc5(z))
        z = SELU(self.fc6(z))
        z = SELU(self.fc7(z))

        # MODIFIED: reshape from 2048 ➝ (32, 16, 2)
        z = z.reshape([batch_size, 32, 16, 2])  # CHANGED

        z = SELU(self.convt0(self.bn4(z)))
        z = self.convt1(self.bn5(z))
        z = self.convt2(self.bn6(z))

        # MODIFIED: reshape output to (128, 16)
        z = z.reshape([batch_size, 128, 16])  # CHANGED
        return z

    def forward(self, x):

        mu, d = self.encoder(x)
        latent_dist = LowRankMultivariateNormal(mu, 0 * mu.unsqueeze(-1), d)
        z = latent_dist.rsample()
        x_rec = self.decoder(z)

        pxz_term = -0.5 * x.shape[1] * (np.log(2 * np.pi / self.model_precision))
        l2s = torch.sum(torch.pow(x.reshape(1, -1) - x_rec.reshape(1, -1), 2), dim=1)
        pxz_term = pxz_term - 0.5 * self.model_precision * torch.sum(l2s)

        elbo = -0.5 * (torch.sum(torch.pow(z, 2)) + self.latent_dim * np.log(2 * np.pi))
        elbo = elbo + pxz_term
        elbo = elbo + torch.sum(latent_dist.entropy())

        return -elbo, mu, d, x_rec
    def train_epoch(self, train_loader):
        self.train()
        train_loss = 0.0
        for _, data in enumerate(train_loader):
            self.optimizer.zero_grad()
            N = data.shape[0]
            data = data.to(self.device)
            loss,_,_,_, = self.forward(data)
            train_loss += loss.item()
            loss.backward()
            self.optimizer.step()
            
        train_loss /= len(train_loader.dataset)
        
        print('Epoch: {} Average loss: {:.4f}'.format(self.epoch, \
                train_loss))
        self.epoch += 1
        return train_loss

    def test_epoch(self, test_loader):

        self.eval()
        test_loss = 0.0
        with torch.no_grad():
            for i, data in enumerate(test_loader):
                
                data = data.to(self.device)
                loss,mu,d,x_rec = self.forward(data)
                test_loss += loss.item()
        
        ind_d=random.randint(0,data.shape[0]-1)
        plt.subplot(1,2,1)
        
        plt.imshow(data.detach().cpu().numpy()[ind_d],origin='lower',cmap = newcmp)
        plt.axis('tight')
        plt.clim([0.3,1])
        plt.subplot(1,2,2)
        plt.imshow(x_rec.detach().cpu().numpy()[ind_d],origin = 'lower',cmap = newcmp)
        plt.axis('tight')
        plt.clim([0.3,1])
        
        display.clear_output(wait=True)
        display.display(plt.gcf())            
        test_loss /= len(test_loader.dataset)
        print('Test loss: {:.4f}'.format(test_loss))
        return test_loss


    def train_loop(self, loaders, epochs=100, test_freq=2, save_freq=10):

        print("="*40)
        print("Training: epochs", self.epoch, "to", self.epoch+epochs-1)
        print("Training set:", len(loaders['train'].dataset))
        print("Test set:", len(loaders['test'].dataset))
        print("="*40)
        # For some number of epochs...
        for epoch in range(self.epoch, self.epoch+epochs):
            # Run through the training data and record a loss.
            loss = self.train_epoch(loaders['train'])
            self.loss['train'][epoch] = loss
            # Run through the test data and record a loss.
            if (test_freq is not None) and (epoch % test_freq == 0):
                loss = self.test_epoch(loaders['test'])
                self.loss['test'][epoch] = loss
            
            

def normalize(spec):
    return (spec-np.min(spec))/(np.max(spec)-np.min(spec))
def get_spec(sig,fs,params):
    '''
    spectrogram
    '''
    nperseg=params['nperseg']
    noverlap=params['noverlap']
    F,T,S = signal.stft(sig,fs,'hann',nperseg,noverlap);
    spec = np.log(np.abs(S) + 1e-12);
    
    if len(T)>1:
        dt=T[2]-T[1];
    else:
        dt = T;
    
    return [spec, dt, F ,T]

def mod_spectrogram(spec_,F,thres,resize_scale):
    ## downsample spectrograms to 128x512
    
    # normalizr amplitude
    spec_ = (spec_-np.min(spec_))/(np.max(spec_)-np.min(spec_))
   
    # only get frequency between 250Hzn and 10000Hz
    
    F_sample = np.linspace(100,10000,resize_scale[0]) 
    T_sample = np.linspace(0,thres,resize_scale[1]) # downsample to 512

    if spec_.shape[1]==thres:
        interp=intp.RegularGridInterpolator([np.linspace(0, spec_.shape[1], spec_.shape[1]), np.linspace(0, F[-1], len(F))], spec_.T)
    else:
        interp_spec= np.zeros((F.shape[0],thres))
        shift = (thres-spec_.shape[1])//2
        
        interp_spec[:,shift:shift+spec_.shape[1]] = spec_
        interp=intp.RegularGridInterpolator([np.linspace(0, interp_spec.shape[1],  interp_spec.shape[1]), np.linspace(0, F[-1], len(F))], interp_spec.T)
    
    xx =np.meshgrid(T_sample, F_sample, indexing='ij')
    sample_points = np.array([xx[0].ravel(), xx[1].ravel()]).T
    
    return interp(sample_points, method='linear').reshape(resize_scale[1],resize_scale[0]).T

#def zero_outside_syllables(spec, T, sylltime, fs):
    mask = np.zeros(spec.shape[1], dtype=bool)  # Mask for time bins

    # Convert each [onset, offset] to time (seconds), then find the matching indices in T
    for onset, offset in sylltime:
        t_start = onset / fs
        t_end = offset / fs
        mask |= (T >= t_start) & (T <= t_end)

    # Apply mask: keep only time bins within any syllable
    spec_filtered = spec.copy()
    spec_filtered[:, ~mask] = 0
    return spec_filtered

def resample_to_new_time(signal, old_time, new_time, kind='linear', fill_value='extrapolate'):
    """
    Resample a signal from old_time to new_time using interpolation.
    
    Parameters:
        signal (ndarray): Shape (..., time), or just (time,)
        old_time (ndarray): Original time points.
        new_time (ndarray): Desired time points.
        kind (str): Type of interpolation: 'linear', 'nearest', 'cubic', etc.
        fill_value: Value to use for extrapolation.
    
    Returns:
        ndarray: Resampled signal with shape (..., len(new_time))
    """
    signal = np.asarray(signal)
    squeezed = False

    if signal.ndim == 1:
        signal = signal[np.newaxis, :]
        squeezed = True

    interp_func = interp1d(old_time, signal, kind=kind, axis=-1, fill_value=fill_value, bounds_error=False)
    resampled = interp_func(new_time)

    if squeezed:
        resampled = resampled[0]

    return resampled

def mov_to_vel_estimate(acc_signal,afs):
    b = firwin(numtaps=481, cutoff=[1, 20], fs=afs, pass_zero=False)  
    vel = filtfilt(b, [1], acc_signal)
    vel = vel-np.mean(vel)
    vel = cumtrapz(vel,tmov,initial=0)
    vel_est = vel-medfilt(vel,kernel_size=2501)
    return vel_est

def ifr_gaussian_kt(filespks, fs, kt, sigma_s, method='auto', support=5, uniform_tol=1e-9):
    """
    Instantaneous firing rate via Gaussian kernel density evaluated on kt.

    Parameters
    ----------
    filespks : array-like or list-of-arrays
        Spike times either as *sample indices* (integers) or *seconds* (floats).
        Integers are interpreted as samples and converted to seconds via fs.
        Floats are treated as seconds.
    fs : float
        Sampling rate (Hz) used if filespks are sample indices.
    kt : array-like
        Time grid in seconds where IFR is evaluated (e.g., 0:binsz/1000:len(S)/fs).
    sigma_s : float
        Gaussian kernel std in seconds (e.g., 0.010 for 10 ms).
    method : {'auto','direct','fft'}
        'direct' = exact KDE sum at each kt; 'fft' = bin+convolve (uniform kt only);
        'auto'   = pick based on size and uniformity.
    support : float
        Kernel half-width in units of sigma for the FFT method (e.g., 5 => ±5σ).
    uniform_tol : float
        Absolute tolerance to treat kt as uniform.

    Returns
    -------
    ifr : ndarray
        Instantaneous firing rate at kt (Hz), same shape as kt.
    """
    def _flat_float(x):
        arr = np.asarray(x, dtype=object)
        if arr.dtype == object:
            parts = []
            for xi in arr.ravel():
                if xi is None:
                    continue
                parts.append(np.asarray(xi).ravel())
            return (np.concatenate(parts).astype(float) if parts else np.array([], dtype=float))
        else:
            return np.asarray(x, dtype=float).ravel()

    # Coerce kt and ensure sorted (monotonic)
    kt_arr = _flat_float(kt)
    if kt_arr.size == 0:
        return np.zeros(0, dtype=float)
    if not np.all(np.diff(kt_arr) >= 0):
        kt_arr = np.sort(kt_arr)

    # Infer dt & uniformity
    if kt_arr.size > 1:
        diffs = np.diff(kt_arr)
        dt = float(np.median(diffs))
        is_uniform = np.allclose(diffs, dt, rtol=0.0, atol=uniform_tol)
    else:
        diffs = np.array([])
        dt = np.nan
        is_uniform = True  # single point

    # Coerce spikes; detect seconds vs samples by integer-ness
    spk = _flat_float(filespks)
    if spk.size == 0:
        return np.zeros_like(kt_arr)

    # If all spikes are (near) integers -> treat as samples; else seconds
    if np.allclose(spk, np.round(spk)):
        spike_t = spk.astype(float) / float(fs)
    else:
        spike_t = spk.astype(float)

    # Choose method
    if method == 'auto':
        cost_direct = kt_arr.size * spike_t.size
        method = 'fft' if (is_uniform and cost_direct > 5e6) else 'direct'

    if (method != 'fft') or (not is_uniform):
        # DIRECT evaluation (exact KDE; units Hz)
        diff = kt_arr[:, None] - spike_t[None, :]
        inv = 1.0 / (sigma_s * np.sqrt(2.0 * np.pi))
        ifr = inv * np.exp(-0.5 * (diff / sigma_s) ** 2).sum(axis=1)
        return ifr

    # FFT / discrete convolution on uniform kt
    start_t = kt_arr[0]
    idx = np.floor((spike_t - start_t) / dt + 0.5).astype(int)
    valid = (idx >= 0) & (idx < kt_arr.size)
    x = np.zeros(kt_arr.size, dtype=float)
    if valid.any():
        np.add.at(x, idx[valid], 1.0)

    half = int(np.ceil(sigma_s * support / dt))
    half = max(1, half)
    tker = (np.arange(-half, half + 1, dtype=float)) * dt
    g = (1.0 / (sigma_s * np.sqrt(2.0 * np.pi))) * np.exp(-0.5 * (tker / sigma_s) ** 2)

    # Discrete conv approximates the integral -> multiply by dt to get Hz
    ifr = dt * np.convolve(x, g, mode='same')
    return ifr


from sklearn.base import clone
from scipy.stats import wilcoxon
import matplotlib.colors as mcolors

def xgb_group_cv_pvalue(
    X, y, groups, model_base,
    n_splits=None,            # if None -> min(max_splits, n_groups)
    max_splits=40,
    make_plot=True,
    ax=None,                  # axis for FULL model plot (optional)
    zero_method="wilcox",
    save_dir=None             # NEW: directory to save figures (optional)
):
    """
    GroupKFold CV; compares FULL model vs train-mean NULL; also fits VOC-only and MOV-only models.
    Returns Wilcoxon p for FULL>NULL and FULL>VOC, FULL>MOD.
    makes scatter plots for VOC-only and MOV-only, and saves figures if save_dir is given.
    """
    groups = np.asarray(groups)
    n_groups = len(np.unique(groups))
    if n_splits is None:
        n_splits = min(max_splits, n_groups)

    gkf = GroupKFold(n_splits=n_splits)

    deltas = []
    deltas_mov = []
    deltas_voc = []
    r2s_mov, r2s_voc, r2s_model = [], [], []
    pct_improvements_full, pct_improvements_voc, pct_improvements_mov = [], [], []

    # splits for your feature subsets (adjust the 32 boundary if needed)
    X_mov_only = X.iloc[:, 32:]
    X_voc_only = X.iloc[:, :32]

    # collect for plots
    y_true_all, y_pred_all, fold_ids = [], [], []
    y_true_all_voc, y_pred_all_voc, fold_ids_voc = [], [], []
    y_true_all_mov, y_pred_all_mov, fold_ids_mov = [], [], []

    for k, (tr, te) in enumerate(gkf.split(X, y, groups)):
        print('Fold '+str(k))
        # safe slicing whether y is np.ndarray or pd.Series
        y_tr = y[tr] if isinstance(y, np.ndarray) else y.iloc[tr]
        y_te = y[te] if isinstance(y, np.ndarray) else y.iloc[te]
        X_tr = X.iloc[tr] if hasattr(X, "iloc") else X[tr]
        X_te = X.iloc[te] if hasattr(X, "iloc") else X[te]

        # null = train-fold mean rate
        y_mean = np.mean(y_tr)
        y_null = np.full_like(y_te, y_mean, dtype=float)

        # ===== FULL model =====
        model = clone(model_base)
        model.fit(X_tr, y_tr)
        y_hat = np.clip(model.predict(X_te), 1e-9, None)

        r2s_model.append(r2_score(y_te, y_hat))
        s_model = -mean_squared_error(y_te, y_hat)
        s_null  = -mean_squared_error(y_te, y_null)
        deltas.append(s_model - s_null)
        pct_improvements_full.append(100.0 * (s_null - s_model) / s_null)

        y_true_all.append(np.asarray(y_te))
        y_pred_all.append(np.asarray(y_hat))
        fold_ids.append(np.full(y_te.shape, k, dtype=int))

        # ===== VOC-only =====
        model_v = clone(model_base)
        model_v.fit(X_voc_only.iloc[tr], y_tr)
        y_hat_v = np.clip(model_v.predict(X_voc_only.iloc[te]), 1e-9, None)

        r2s_voc.append(r2_score(y_te, y_hat_v))
        # Δ(full over voc-only): (-MSE_full) - (-MSE_voc) = MSE_voc - MSE_full
        deltas_voc.append(s_model - (-mean_squared_error(y_te, y_hat_v)))
        pct_improvements_voc.append(100.0 * (s_null - (-mean_squared_error(y_te, y_hat_v))) / s_null)

        y_true_all_voc.append(np.asarray(y_te))
        y_pred_all_voc.append(np.asarray(y_hat_v))
        fold_ids_voc.append(np.full(y_te.shape, k, dtype=int))

        # ===== MOV-only =====
        model_m = clone(model_base)
        model_m.fit(X_mov_only.iloc[tr], y_tr)
        y_hat_m = np.clip(model_m.predict(X_mov_only.iloc[te]), 1e-9, None)

        r2s_mov.append(r2_score(y_te, y_hat_m))
        deltas_mov.append(s_model - (-mean_squared_error(y_te, y_hat_m)))
        pct_improvements_mov.append(100.0 * (s_null - (-mean_squared_error(y_te, y_hat_m))) / s_null)

        y_true_all_mov.append(np.asarray(y_te))
        y_pred_all_mov.append(np.asarray(y_hat_m))
        fold_ids_mov.append(np.full(y_te.shape, k, dtype=int))

    deltas     = np.asarray(deltas)
    deltas_mov = np.asarray(deltas_mov)
    deltas_voc = np.asarray(deltas_voc)

    # Wilcoxon tests
    stat, p = wilcoxon(deltas, alternative="greater", zero_method=zero_method)          # FULL > NULL
    _, p_over_voc = wilcoxon(deltas_voc, alternative="greater")                         # FULL > VOC-only
    _, p_over_mov = wilcoxon(deltas_mov, alternative="greater")                         # FULL > MOV-only

    # =======================
    # Plot helpers
    # =======================
    def _scatter_pred_vs_actual_row(
        y_trues_list, y_preds_list, folds_list,
        titles, n_splits, save_path=None
    ):
        fig = plt.figure(figsize=(18, 6))
        gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 0.05])  # last slot = cbar
        axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
        cax = fig.add_subplot(gs[0, 3])  # explicit colorbar axis
    
        norm = mcolors.Normalize(vmin=0, vmax=n_splits - 1)
        cmap = cm.get_cmap("tab20b", n_splits)
    
        # global min/max for consistent axes
        all_true = np.concatenate(sum(y_trues_list, []))
        all_pred = np.concatenate(sum(y_preds_list, []))
        lo = float(np.min([all_true.min(), all_pred.min()]))
        hi = float(np.max([all_true.max(), all_pred.max()]))
    
        for ax, y_trues, y_preds, folds, title in zip(axes, y_trues_list, y_preds_list, folds_list, titles):
            y_true_concat = np.concatenate(y_trues)
            y_pred_concat = np.concatenate(y_preds)
            fold_concat   = np.concatenate(folds)
    
            for k in range(n_splits):
                m = (fold_concat == k)
                if np.any(m):
                    ax.scatter(y_true_concat[m], y_pred_concat[m],
                               s=5, alpha=0.5, color=cmap(norm(k)))
            ax.plot([lo, hi], [lo, hi], linestyle=":", linewidth=3, color="red")
            ax.set_title(title)
            ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal", "box")
            ax.set_xlabel("Actual FR (Hz)")
    
        # only left plot has y-label
        axes[0].set_ylabel("Predicted FR (Hz)")
    
        # shared colorbar on far right
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        fig.colorbar(sm, cax=cax, label="Fold index")
    
        fig.tight_layout()
        if save_path is not None:
            fig.savefig(save_path, dpi=200, bbox_inches="tight")
        return fig, axes


    # =======================
    # Make/summarize plots
    # =======================
    fig_full = None
    fig_voc  = None
    fig_mov  = None
    saved_full = saved_voc = saved_mov = None

    if make_plot:
        titles = [
            f"Full (p={p:.3g})",
            f"Voc-only (p={p_over_voc:.3g})",
            f"Mov-only (p={p_over_mov:.3g})"
        ]
        fig_row, axes_row = _scatter_pred_vs_actual_row(
            [y_true_all, y_true_all_voc, y_true_all_mov],
            [y_pred_all, y_pred_all_voc, y_pred_all_mov],
            [fold_ids, fold_ids_voc, fold_ids_mov],
            titles, n_splits,
            save_path=(Path(save_dir)/"scatter_row.png" if save_dir else None)
        )
    fig_boxes = None
    saved_boxes = None
    
    # now the box plots
    if make_plot:
        labels = ["Full", "Voc-only", "Mov-only"]
        all_r2_values = [np.asarray(r2s_model), np.asarray(r2s_voc), np.asarray(r2s_mov)]
        all_pct_values = [
            np.asarray(pct_improvements_full),
            np.asarray(pct_improvements_voc),
            np.asarray(pct_improvements_mov),
        ]
    
        fig_boxes, axes_boxes = plt.subplots(1, 2, figsize=(10, 6))
        fig_boxes.subplots_adjust(wspace=0.5)
        # Shared style dict 
        style = dict(
            patch_artist=True,
            boxprops=dict(facecolor='lightblue', linewidth=2),
            whiskerprops=dict(linewidth=2),
            capprops=dict(linewidth=2),
            medianprops=dict(linewidth=2),
            flierprops=dict(marker='o', color='red', markersize=6, linewidth=2)
        )
    
        # Left: R^2 boxplot
        axes_boxes[0].boxplot(all_r2_values, **style)
        axes_boxes[0].set_xticks([1, 2, 3], labels=labels)
        axes_boxes[0].set_ylabel(f"Mean R\u00b2")
        axes_boxes[0].set_ylim(-1,1)

    
        # Right: % improvement vs null (higher is better)
        axes_boxes[1].boxplot(all_pct_values, **style)
        axes_boxes[1].set_xticks([1, 2, 3], labels=labels)
        axes_boxes[1].set_ylabel("% improvement from null")


    # Save
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        saved_boxes = os.path.join(save_dir, "boxplots_row.png")
        fig_boxes.savefig(saved_boxes, dpi=200, bbox_inches="tight")


    return {
        "deltas": deltas,
        "pct_imp_full": pct_improvements_full,
        "pct_imp_voc": pct_improvements_voc,
        "pct_imp_mov": pct_improvements_mov,
        "r2s_model": np.asarray(r2s_model),
        "r2s_voc": np.asarray(r2s_voc),
        "r2s_mov": np.asarray(r2s_mov),
        "frac_groups_positive": float(np.mean(deltas > 0)),
        "wilcoxon_stat": float(wilcoxon(deltas, alternative='greater', zero_method=zero_method)[0]),
        "p_value": float(p),                 # FULL > NULL
        "p_over_voc": float(p_over_voc),    # FULL > VOC-only
        "p_over_mov": float(p_over_mov),    # FULL > MOV-only
        "n_splits": int(n_splits),
        "n_groups": int(n_groups),
        "fig_full": fig_full, "fig_voc": fig_voc, "fig_mov": fig_mov,
        "saved_full": str(saved_full) if save_dir is not None else None,
        "saved_voc":  str(saved_voc)  if save_dir is not None else None,
        "saved_mov":  str(saved_mov)  if save_dir is not None else None,
    }

#####################################
params = {
    'nperseg': 128,
    'noverlap': 128/2,
}

thres = 16

# === Custom Jet Colormap (First 20 colors black) ===
jet_colors = cm.get_cmap('jet', 256)(np.linspace(0, 1, 256))
jet_colors[:20] = np.array([0, 0, 0, 1])  # Set first 20 entries to black
newcmp = ListedColormap(jet_colors)

# === Suppress Warnings ===
warnings.filterwarnings("ignore")

sampling_rate = 20000
fir_coeff = signal.firwin(numtaps=101,cutoff=[150,9999],fs=sampling_rate,pass_zero=False)

#================================================ vibin' in this bih

# save path and dbase list 
base_savepath = os.path.join('X:\\Budgie','AAC_sorted_dbases_cbj','full_segments','infoadded','NEW_whisperseg','joint_modeling_ifr_final')

dbase_list = glob.glob(os.path.join('X:\\Budgie','AAC_sorted_dbases_cbj','full_segments','infoadded','NEW_whisperseg','*.mat'))

for yeet,pick_dbase in enumerate(dbase_list[84:],start=84):
    print('dbase '+str(yeet))
    #pick_dbase = 'dbase0880_20190831_segs_vocsegs_chan16_FIR_zz_segs.mat'
    #full_path = os.path.join('X:\\Budgie','AAC_sorted_dbases_cbj','full_segments','infoadded','NEW_whisperseg',pick_dbase)
    data = sio.loadmat(pick_dbase,struct_as_record=False,squeeze_me=True)
    dbase = data['dbase']
    
    # get path to dbase txt files, 
    sortfiles = dbase.sortedfiles
    thispath = dbase.PathName
    if not hasattr(dbase,'headbob'):
        hbfiles = []
    else:  
        hbfiles = [i for i, arr in enumerate(dbase.headbob.boutTimes) if arr.size >0]
    hbfiles = np.array(hbfiles)+1 #change back to matlab index
    shbfiles = np.intersect1d(hbfiles,sortfiles)
    if not any(shbfiles):
        continue
    num_pairs = sum(arr.shape[0] if arr.ndim == 2 else 1 for arr in dbase.headbob.boutTimes[shbfiles-1] if arr.size)
    if num_pairs<20:
        continue
    
    # load bird's model
    birdID = dbase.birdID
    base_savename = birdID+'_'+dbase.date+'_'+dbase.ephyschan
    this_savepath = os.path.join(base_savepath,base_savename)
    os.makedirs(this_savepath, exist_ok=True)
    model_path = os.path.join('X:\\Budgie','AAC_sorted_dbases_cbj','full_segments','infoadded','NEW_whisperseg','VAEs','VAE_for_hb'+birdID+'_bin_'+str(thres)+'.pth')
    model= VAE(latent_dim=32,model_precision=10,lr =1e-4)
    model.load_state_dict(torch.load(model_path, map_location=torch.device('cpu')))
    model.eval()
    
    hb_pad = 0#.5 # peri-hb pad to include in s 
    afs = 5000
    

    all_latents = []
    all_spks_25 = []
    all_spks_50 = []
    all_spks_100 = []
    all_spks_200 = []
    all_ifrs = []
    all_mov = []
    all_z_vel = []
    all_x_vel = []
    all_y_vel = []
    all_hb_time = []
    all_abs_hb_time = []
    all_hb_phase1, all_hb_phase2 = [],[]
    all_hb_time_mins = []
    all_is_silence = []
    all_file_ID = []
    
    # averaged move features instead
    all_m_z_acc = []
    all_m_y_acc = []
    all_m_x_acc = []
    all_m_x_vel = []
    all_m_z_vel_1,all_m_z_vel_2,all_m_z_vel_3,all_m_z_vel_4 = [],[],[],[]
    all_m_x_vel_1,all_m_x_vel_2,all_m_x_vel_3,all_m_x_vel_4 = [],[],[],[]
    all_m_y_vel = []
    # all_m_z_hil_i,all_m_z_hil_2,all_m_z_hil_3,all_m_z_hil_4 = [],[],[],[]
    # all_m_y_hil = []
    # all_m_x_hil_i,all_m_x_hil_2,all_m_x_hil_3,all_m_x_hil_4 = [],[],[],[]
    
    
    
    
    for i,num in enumerate(shbfiles):
        print('File '+str(i+1)+' of '+str(len(shbfiles)))
        sfilename = dbase.SoundFiles[num-1].name
        zfilename = sfilename[:sfilename.find('0.txt')]+'18.txt'
        xfilename = sfilename[:sfilename.find('0.txt')]+'19.txt'
        yfilename = sfilename[:sfilename.find('0.txt')]+'20.txt'
        
        # spikes
        if not hasattr(dbase.sortedspks[dbase.sortedfiles==num][0],'shape'):
            spks = dbase.sortedspks[dbase.sortedfiles==num]
            sspks = dbase.sortedspks[dbase.sortedfiles==num]
        else:
            spks = dbase.sortedspks[dbase.sortedfiles==num][0]
            sspks = dbase.sortedspks[dbase.sortedfiles==num][0]
        spks = spks/sampling_rate
        # headbob
    
        bobtimes = dbase.headbob.indbobTimes[num-1]
        bouttimes = dbase.headbob.boutTimes[num-1]
        if bouttimes.ndim==1:
            bouttimes= bouttimes.reshape(-1,2)/sampling_rate
        else:
            bouttimes = bouttimes/sampling_rate
        if bobtimes.ndim==1:
            bobtimes = bobtimes.reshape(-1,2)/sampling_rate
        else:
            bobtimes = bobtimes/sampling_rate
            
    
        if any((bobtimes[:,1]-bobtimes[:,0])<0.05):
            print("WARNING BOB TIMES MAYBE MESSED UP, File "+str(num))
        #sound
        sound = np.loadtxt(os.path.join(thispath,sfilename),skiprows=3)
        sound = signal.filtfilt(fir_coeff,1,sound)
        spec, dt, F, T = get_spec(sound, sampling_rate, params)
        spec = normalize(spec)
        # if dbase.voc_segs_zz.whispersegs[dbase.voc_segs_zz.ifile==num].size==0:
        #     continue
        # else:
        #     sylltime = dbase.voc_segs_zz.whispersegs[dbase.voc_segs_zz.ifile== num][0]
        # if len(sylltime)<3:
        #     continue
        # spec_segtimes= np.round(sylltime / sampling_rate / dt).astype(int)
        # labs = dbase.whisper_seg_labels_zz[dbase.voc_segs_zz.ifile== num][0]
        #warbpass = dbase.voc_segs_zz.warb.warbpass[dbase.voc_segs_zz.ifile== num][0]
        
        # load in acceleromter channels
        movz = np.loadtxt(os.path.join(thispath,zfilename),skiprows=3)
        movx = np.loadtxt(os.path.join(thispath,xfilename),skiprows=3)
        #movy = np.loadtxt(os.path.join(thispath,yfilename),skiprows=3)
        tmov = np.linspace(0, len(movz)/5000, len(movz))
        
        # movz_filt = medfilt(movz,kernel_size=501)-medfilt(movz,kernel_size=2501)
        # movz_filt = savgol_filter(movz_filt,window_length=101,polyorder=2)
        # ds_mov = resample_to_new_time(movz_filt, tmov, T)
        
        # velocity estimates
        vel_z = mov_to_vel_estimate(movz,afs)
        ds_velz_est = resample_to_new_time(vel_z, tmov, T)
        #vel_y = mov_to_vel_estimate(movy, afs)
        #ds_vely_est = resample_to_new_time(vel_y, tmov, T)
        vel_x = mov_to_vel_estimate(movx, afs)
        ds_velx_est = resample_to_new_time(vel_x, tmov, T)
    
        #ifr
        # ifr = []
        # ifr = ifr_gaussian_kt(sspks,sampling_rate,T,0.1)
        
        
        # assign time in hb and phase for each point in T
        phases = np.full_like(T,np.nan)
        for onset,offset in bobtimes: # assign time in hb phase
            in_cyc = (T>=onset)&(T<=offset)
            t_in = T[in_cyc]
            norm_pos = (t_in-onset)/(offset-onset)
            phases[in_cyc] = norm_pos*2*np.pi
            
        hb_times = np.full_like(T,np.nan)          
        for onset,offset in bouttimes: # assign time in hb
            in_hb = (T>=onset)&(T<=offset)
            t_in = T[in_hb]                    
            hb_times[in_hb] = (t_in-onset)/(offset-onset)
         
        # populate the latents for each headbob
        for k,(onset,offset) in enumerate(bouttimes):
            relspks = spks-onset-hb_pad
            T_inds = np.where((T>=onset-hb_pad)&(T<=offset+hb_pad))[0]
            len_T = len(T_inds)
            samples = int(max(0,(len_T-thres-1)))
            this_spec = spec[:,T_inds]
            #this_mov = ds_mov[T_inds]
            this_z_vel = np.round(ds_velz_est[T_inds]*1000,3)
            this_z_hil = np.abs(hilbert(this_z_vel))
            #this_y_vel = np.round(ds_vely_est[T_inds]*1000,3)
            this_x_vel = np.round(ds_velx_est[T_inds]*1000,3)
            this_x_hil = np.abs(hilbert(this_z_vel))
            these_hb_times= hb_times[T_inds]
            these_abs_times = T[T_inds]-onset
            these_phases = phases[T_inds]
            
            prev15 = np.arange(T_inds[0]-15,T_inds[0],dtype=T_inds.dtype)
            ifr_T = np.concatenate([prev15,T_inds])
            these_ifr = ifr_gaussian_kt(sspks, sampling_rate, T[ifr_T], 0.01)
            
            if samples>0:
    
                sample_reg = np.arange(0,len_T-thres-1)
                for j in range(len(sample_reg)):
                    idx_tmp = sample_reg[j]
                    these_relspks = relspks-idx_tmp*dt
                    this_chunk = this_spec[:,idx_tmp:min(len_T,idx_tmp+thres)]
                    this_chunk = mod_spectrogram(this_chunk, F, thres, (128,thres))
                    this_chunk = np.asarray(this_chunk, dtype='float32')
                    spec_tens = torch.tensor(this_chunk,dtype=torch.float32).unsqueeze(0)
                    with torch.no_grad():
                        _,mu,d,recon = model.forward(spec_tens)
                    all_latents.append(mu.numpy()[0])
                    all_spks_25.append(np.sum((these_relspks>-0.025)&(these_relspks<-0.005)))
                    all_spks_50.append(np.sum((these_relspks>-0.05)&(these_relspks<-0.005)))
                    all_spks_100.append(np.sum((these_relspks>-0.1)&(these_relspks<-0.005)))
                    all_spks_200.append(np.sum((these_relspks>-0.2)&(these_relspks<-0.005)))
                    if any(ifr_T[idx_tmp:idx_tmp+15])<0:
                        continue
                    all_ifrs.append(np.mean(these_ifr[idx_tmp:idx_tmp+15]))#
                    #all_mov.append(this_mov[idx_tmp:idx_tmp+thres])
                    # all_z_vel.append(this_z_vel[idx_tmp:idx_tmp+thres])
                    # all_y_vel.append(this_y_vel[idx_tmp:idx_tmp+thres])
                    # all_x_vel.append(this_x_vel[idx_tmp:idx_tmp+thres])
                    all_hb_phase1.append(these_phases[idx_tmp])
                    all_hb_phase2.append(these_phases[idx_tmp+thres])
                    all_hb_time.append(np.mean(these_hb_times[idx_tmp:idx_tmp+thres]))
                    all_abs_hb_time.append(np.mean(these_abs_times[idx_tmp:idx_tmp+thres]))
                    
                    
                    # do averaged features on a half thres basis, getting the before, during, and after thres/2 bins
                    
                    # this bin first 
                    #all_m_z_acc.append(np.mean(this_mov[idx_tmp:idx_tmp+thres]))
                    all_m_z_vel_1.append(np.mean(this_z_vel[idx_tmp]))
                    all_m_x_vel_1.append(np.mean(this_x_vel[idx_tmp]))
                    # all_m_z_hil_i.append(np.mean(this_z_hil[idx_tmp:idx_tmp+thres/2]))
                    # all_m_x_hil_i.append(np.mean(this_x_hil[idx_tmp:idx_tmp+thres/2]))
                    # this bin last
                    all_m_z_vel_2.append(np.mean(this_z_vel[idx_tmp+thres]))
                    all_m_x_vel_2.append(np.mean(this_x_vel[idx_tmp+thres]))
                    # all_m_z_hil_2.append(np.mean(this_z_hil[idx_tmp+thres/2:idx_tmp+thres]))
                    # all_m_x_hil_2.append(np.mean(this_x_hil[idx_tmp+thres/2:idx_tmp+thres]))
                    
                    # prev bin
                    # all_m_z_vel_3.append(np.mean(this_z_vel[idx_tmp-thres/2:idx_tmp]))
                    # all_m_x_vel_3.append(np.mean(this_x_vel[idx_tmp-thres/2:idx_tmp]))
                    # all_m_z_hil_3.append(np.mean(this_z_hil[idx_tmp-thres/2:idx_tmp]))
                    # all_m_x_hil_3.append(np.mean(this_x_hil[idx_tmp-thres/2:idx_tmp]))
                    
                    # prev prev bin
                    # all_m_z_vel_4.append(np.mean(this_z_vel[idx_tmp-thres:idx_tmp-thres/2]))
                    # all_m_x_vel_4.append(np.mean(this_x_vel[idx_tmp-thres:idx_tmp-thres/2]))
                    # all_m_z_hil_4.append(np.mean(this_z_hil[idx_tmp-thres:idx_tmp-thres/2]))
                    # all_m_x_hil_4.append(np.mean(this_x_hil[idx_tmp-thres:idx_tmp-thres/2]))
                    
                    #all_m_x_acc.append(np.mean(this_x_acc[idx_tmp:idx_tmp+thres]))
                    #all_m_y_acc.append(np.mean(this_y_acc[idx_tmp:idx_tmp+thres]))
                    
                    all_file_ID.append(i)
                    
                
    
    ### ======= now we can build that mf model
    if sum(all_spks_50)<20:
        continue
    import pandas as pd
    import xgboost as xgb
    from sklearn.metrics import r2_score,make_scorer,mean_squared_error
    from sklearn.model_selection import GroupKFold
    from scipy.stats import wilcoxon
    # Set global font size
    matplotlib.rc('font', size=16)
    
    X_array = np.vstack(all_latents)
    # column names
    latent_cols = [f'latent_{i}' for i in range(all_latents[0].shape[0])]
    # = [f'mov_{i}' for i in range(4)]# TIMES 3 if using all three axes all_z_vel[0].shape[0]*3
    all_cols = latent_cols# + mov_cols
    
    #create pandas dataframe
    X_all = pd.DataFrame(X_array, columns = all_cols)
    X_all['hb_phase_1'] = all_hb_phase1
    X_all['hb_phase_2'] = all_hb_phase2
    X_all['hb_time'] = all_hb_time
    X_all['hb_abs_time'] = all_abs_hb_time
    #X_all['m_z_acc'] = all_m_z_acc
    X_all['m_z_vel_1'] = all_m_z_vel_1
    X_all['m_z_vel_2'] = all_m_z_vel_2
    # X_all['m_z_vel_3'] = all_m_z_vel_3
    # X_all['m_z_vel_4'] = all_m_z_vel_4
    # X_all['m_z_hil_i'] = all_m_z_hil_i
    # X_all['m_z_hil_2'] = all_m_z_hil_2
    # X_all['m_z_hil_3'] = all_m_z_hil_3
    # X_all['m_z_hil_4'] = all_m_z_hil_4
    X_all['m_x_vel_1'] = all_m_x_vel_1
    X_all['m_x_vel_2'] = all_m_x_vel_2
    # X_all['m_x_hil_i'] = all_m_x_hil_i
    # X_all['m_x_hil_2'] = all_m_x_hil_2
    # X_all['m_x_hil_3'] = all_m_x_hil_3
    # X_all['m_x_hil_4'] = all_m_x_hil_4
    y_all = np.array(all_ifrs) # PICK spike bin =====
    
    #train test split by file ID
    groups = np.array(all_file_ID)
    
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
    
    # Run through all validation and model splits 
    out = xgb_group_cv_pvalue(X_all,y_all,groups,model_base,save_dir = this_savepath)

    # save data 
    with open(os.path.join(this_savepath, 'summary_dat.pkl'), "wb") as f:
        pickle.dump(out, f)
    
   
    
