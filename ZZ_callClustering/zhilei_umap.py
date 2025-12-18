import shutil, os, sys, pickle, tqdm, math, netCDF4
from scipy.io import wavfile, loadmat, netcdf_file
import scipy
import pandas as pd
import numpy as np
import h5py
import torch
from torch.utils.data import Dataset, DataLoader
import matplotlib.pyplot as plt
import librosa
from joblib import Parallel, delayed
from scipy.interpolate import RegularGridInterpolator
from scipy.signal import butter, lfilter
import warnings


def ZZ_get_spec(s1, s2, audio, fs, nperseg=256, noverlap=128, target_freqs=None, target_time_dim=None, bool_clip=True, clim=(-13,-5)):
    # function modified from Goffinet/Pearson elife paper to make spectrograms
    # bool_clip: whether to normalize and clip spectrograms
    # calculate spectrogram on the audio segment
    epsilon = 10**clim[0]  # add a small offset when calculating spectrogram to avoid np.log(0)
    temp_audio = audio[max(0,s1):min(len(audio),s2)]
    f, t, spec = scipy.signal.stft(temp_audio, fs=fs, nperseg=nperseg, noverlap=noverlap)
    spec = np.log(np.abs(spec) + epsilon)
    # interpolate into desired dimension if specified 
    if target_freqs is None:
        target_freqs = f
    if target_time_dim is None:
        target_times = t
    else:
        target_times = np.linspace(0, t[-1], target_time_dim)
    fv, tv = np.meshgrid(target_freqs, target_times, indexing='ij')
    interp = scipy.interpolate.RegularGridInterpolator((f, t), spec, bounds_error=False, fill_value=clim[0])
    points = np.array([fv.ravel(), tv.ravel()]).T
    interp_spec = interp(points)
    interp_spec = interp_spec.reshape(fv.shape)
    # normalize and clip if specified
    if bool_clip:
        interp_spec -= clim[0]
        interp_spec /= (clim[1]-clim[0])
        interp_spec = np.clip(interp_spec, 0.0, 1.0)
    return interp_spec, t[1], target_freqs, target_times


def ZZ_specFromWavZZ_v1(fn, p, syl):
    # calculate spectrograms from wav file
    # check if annotation file exist
    fn_label = fn.replace('.wav','.label.txt')
    # given wav file and selected syllables, output spectrograms and meta info
    specs = []
    info = pd.DataFrame()
    # when calcualte spectrogram, include a small padding to avoid artifacts at the boundary
    if os.path.exists(fn_label) and os.path.getsize(fn_label)>0:
        labels = np.genfromtxt(fn_label, dtype=str, delimiter=None, encoding='utf-8-sig')
        labels = np.atleast_1d(labels)  # deal with only one segment in a file
        # check if it has target syllable
        idx = [ii for ii in range(labels.shape[0]) if labels[ii] in syl]
        if len(idx)>0:
            fn_time = fn_label.replace('.label.txt', '.time.txt')
            seg = np.loadtxt(fn_time, delimiter=',', dtype=int)
            seg = np.atleast_2d(seg)
            # load wav file
            audio, fs = librosa.load(fn, sr=None)
            # butterworth filter
            audio = butter_bandpass_filter(audio, p['target_freqs'][0], p['target_freqs'][-1], fs, order=5)
            for m in idx:
                # calculate spectrogram
                spec, dt, target_freqs, target_times = ZZ_get_spec(seg[m,0], seg[m,1], audio, fs, nperseg=p['nperseg'], noverlap=p['noverlap'], 
                                                       target_freqs=p['target_freqs'], target_time_dim=p['target_time_dim'],
                                                       bool_clip=p['bool_clip'], clim=p['clim'])
                specs.append(spec)
                # record the meta info of this syllable
                row = pd.DataFrame([{'fn_wav':fn, 's_idx':m, 'istart':seg[m,0], 'iend':seg[m,1], 'label':labels[m], 'spec_f':target_freqs, 'spec_t':target_times}])
                info = pd.concat([info, row], ignore_index=True)
    return(specs, info)


def butter_bandpass(lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype="band")
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=3):
    if highcut > int(fs / 2):
        warnings.warn("Highcut is too high for bandpass filter. Setting to nyquist")
        highcut = int(fs / 2)
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def pad_spectrogramZZ(spectrogram, pad_length):
    """ Pads a spectrogram to being a certain length, if spectrogram is long, cutoff to the center part
    """
    if spectrogram.shape[1]<=pad_length:
        excess_needed = pad_length - np.shape(spectrogram)[1]
        pad_left = np.floor(float(excess_needed) / 2).astype("int")
        pad_right = np.ceil(float(excess_needed) / 2).astype("int")
        return np.pad(
            spectrogram, [(0, 0), (pad_left, pad_right)], "constant", constant_values=0
        )
    else:
        to_remove = int((spectrogram.shape[1] - pad_length)/2)
        spec_p = np.zeros([spectrogram.shape[0], pad_length], dtype=spectrogram.dtype)
        spec_p = spectrogram[:, to_remove:(to_remove+pad_length)]
        return(spec_p)