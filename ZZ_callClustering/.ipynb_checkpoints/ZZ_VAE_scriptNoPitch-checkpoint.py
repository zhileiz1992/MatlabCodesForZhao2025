## use Tim Sainburg's segmentation code
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

from vocalseg.dynamic_thresholding import dynamic_threshold_segmentation
from vocalseg.utils import butter_bandpass_filter, spectrogram, int16tofloat32, plot_spec

VAE_script_path = '/home/zz367/ProjectsU/WarbleAnalysis/Jupyter/ZZ_VAE_UMAP/'
sys.path.insert(1, VAE_script_path)
# load Han's script on VAE architecture
import Han_VAE_scriptNoPitch
import Han_VAE_modelNew

num_cpu = 48   # number of cpu threads to use for parallel processing


def ZZ_timSegDynamic_warble(fn_pair):
    # fn_pair = [fn_wav, fd_target, focal_ch, param_seg]
    # first save focal_ch of wav file fn_wav to the target folder fd_target
    # read in wav file fn_wav, segment syllables in channel focal_ch from silence
    # write the results into the fd_target folder: original wav, time.txt, label.txt
    [fn_wav, fd_target, focal_ch, param_seg] = fn_pair
    rate, data = wavfile.read(fn_wav)
    # pick the specific channel 
    if len(data.shape)>1:
        data = data[:, focal_ch]
    # save focal_ch as a separate wav files
    fn_save = os.path.join(fd_target, os.path.basename(fn_wav))
    wavfile.write(fn_save, rate, data)
    # change to float type
    data = data / 32767
    # band pass raw data
    data = butter_bandpass_filter(data, param_seg['spectral_range'][0], param_seg['spectral_range'][1], rate, order=2)
    # segment using the dynamic thresholding 
    results = dynamic_threshold_segmentation(data, rate, **param_seg)
    # save segment times
    pd_time = pd.DataFrame()
    pd_time['onset'] = [int(aa*rate) for aa in results['onsets']]
    pd_time['offset'] = [int(aa*rate) for aa in results['offsets']]
    fn_time = os.path.join(fd_target, os.path.basename(fn_wav).replace('wav', 'time.txt'))
    pd_time.to_csv(fn_time, header=None, index=False)
    # save pseudo-labels
    fn_label = os.path.join(fd_target, os.path.basename(fn_wav).replace('wav', 'label.txt'))
    label_str = ','.join(['a'] * len(pd_time['onset']))
    with open(fn_label, 'w') as f:
        f.write(label_str)

def ZZ_timSegDynamic_warble_v2(fn_pair):
    # fn_pair = [fn_wav, fd_target, focal_ch, param_seg]
    # first save focal_ch of wav file fn_wav to the target folder fd_target
    # read in wav file fn_wav, segment syllables in channel focal_ch from silence
    # write the results into the fd_target folder: original wav, time.txt, label.txt
    [fn_wav, fd_target, focal_ch, param_seg] = fn_pair
    rate, data = wavfile.read(fn_wav)
    # pick the specific channel 
    if len(data.shape)>1:
        data = data[:, focal_ch]
    # save focal_ch as a separate wav files
    fn_save = os.path.join(fd_target, os.path.basename(fn_wav))
    wavfile.write(fn_save, rate, data)
    # change to float type
    data = data / 32767
    # band pass raw data
    data = butter_bandpass_filter(data, param_seg['spectral_range'][0], param_seg['spectral_range'][1], rate, order=2)
    # segment using the dynamic thresholding 
    results = dynamic_threshold_segmentation(data, rate, **param_seg)
    if results:
        # save segment times
        pd_time = pd.DataFrame()
        pd_time['onset'] = [int(aa*rate) for aa in results['onsets']]
        pd_time['offset'] = [int(aa*rate) for aa in results['offsets']]
        fn_time = os.path.join(fd_target, os.path.basename(fn_wav).replace('wav', 'time.txt'))
        pd_time.to_csv(fn_time, header=None, index=False)
        # save pseudo-labels
        fn_label = os.path.join(fd_target, os.path.basename(fn_wav).replace('wav', 'label.txt'))
        label_str = ','.join(['a'] * len(pd_time['onset']))
        with open(fn_label, 'w') as f:
            f.write(label_str)


def ZZ_chopSpectWins_v1(fn_pair):
    # chop spectrograms into windows
    # fn_pair=[fn_wav, focal_ch, param_chop]
    # assume that XXX.time.txt and XXX.label.txt exist for fn_wav in the same folder
    # save results as h5 files for future reuse
    [fn_wav, focal_ch, param_chop] = fn_pair
    fd_save_temp = fn_wav.replace('.wav', '')
    # clear previous results
    if os.path.exists(fd_save_temp):
        shutil.rmtree(fd_save_temp)
    os.makedirs(fd_save_temp)
    delta = param_chop['delta']  # this is the width of spectrogram window, in unit of spectrogram frames, Han uses 32
    hop = param_chop['hop']  # default hop length is 1 frame, convert to sec: hop * dt, Han uses 1
    param_spec = {'nperseg': param_chop['nperseg'], 'noverlap':param_chop['noverlap']}
    # load wav file
    fs, data = wavfile.read(fn_wav)
    if len(data.shape)>1:
        data = data[:, focal_ch]
    data = data / 32767  # convert to float
    # band pass raw data
    data = butter_bandpass_filter(data, param_chop['bp_low'], param_chop['bp_high'], fs, order=2)
    # get spectrograms using Han's method
    [spec_full, dt, F ,T] = Han_VAE_scriptNoPitch.get_spec(data, fs, param_spec)
    # load the segmentation
    seg = np.loadtxt(fn_wav.replace('wav', 'time.txt'), delimiter=',').astype('int')  # load the segmentation
    seg = np.round(seg/fs/dt).astype('int')  # convert syllable boundary to spectrogram frames
    seg[seg>=len(T)]= len(T)
    # load the syllable labels
    label = np.loadtxt(fn_wav.replace('wav', 'label.txt'), dtype='str',delimiter=',')
    # loop through syllables, chop into spectrogram windows
    for n in range(seg.shape[0]):
        spec_tmp = spec_full[:,seg[n,0]:seg[n,1]]
        spec_wins = []
        for t in range(0, spec_tmp.shape[1]-delta, hop):
            interp_spec = scipy.interpolate.RectBivariateSpline(F[:128],np.linspace(0,T[delta],delta), spec_tmp[:128,t:t+delta])
            spec_intp = interp_spec(F[:128],np.linspace(0,T[delta],delta))
            # flip the freq axis so lower freq is at bottom
            spec_intp = np.flip(spec_intp, axis=0)
            # spec_wins.append(spec_intp)
            # save as spectrogram window as h5 file: fn = wavname_syllableNum_windowNum
            fn_h5 = os.path.join(fd_save_temp, f'{os.path.basename(fd_save_temp)}_{label[n]}_{n:05}_{t:05}.h5')
            with h5py.File(fn_h5, 'w') as hf:
                hf.create_dataset('spec', data=spec_intp)
                hf.create_dataset('seg', data=seg[n,:])
                hf.create_dataset('seg_sec', data=seg[n,:]*dt)


def ZZ_chopSpectWins_v2(fn_pair):
    # chop spectrograms into windows
    # fn_pair=[fn_wav, param_chop]
    # assume that XXX.time.txt and XXX.label.txt exist for fn_wav in the same folder
    # return a list of spectrogram windows and associated ids for each window
    [fn_wav, param_chop] = fn_pair
    cond_name = fn_wav.split('/')[-2]
    delta = param_chop['delta']  # this is the width of spectrogram window, in unit of spectrogram frames, Han uses 32
    hop = param_chop['hop']  # default hop length is 1 frame, convert to sec: hop * dt, Han uses 1
    param_spec = {'nperseg': param_chop['nperseg'], 'noverlap':param_chop['noverlap']}
    # load wav file
    fs, data = wavfile.read(fn_wav)
    data = data / 32767  # convert to float
    # band pass raw data
    data = butter_bandpass_filter(data, param_chop['bp_low'], param_chop['bp_high'], fs, order=2)
    # get spectrograms using Han's method
    [spec_full, dt, F ,T] = Han_VAE_scriptNoPitch.get_spec(data, fs, param_spec)
    # load the segmentation
    seg = np.loadtxt(fn_wav.replace('wav', 'time.txt'), delimiter=',').astype('int')  # load the segmentation
    seg = np.round(seg/fs/dt).astype('int')  # convert syllable boundary to spectrogram frames
    seg[seg>=len(T)]= len(T)
    # load the syllable labels
    label = np.loadtxt(fn_wav.replace('wav', 'label.txt'), dtype='str',delimiter=',')
    # loop through syllables, chop into spectrogram windows
    spec_wins = []
    id_wins = []
    for n in range(seg.shape[0]):
        spec_tmp = spec_full[:,seg[n,0]:seg[n,1]]
        for t in range(0, spec_tmp.shape[1]-delta, hop):
            interp_spec = scipy.interpolate.RectBivariateSpline(F[:128],np.linspace(0,T[delta],delta), spec_tmp[:128,t:t+delta])
            spec_intp = interp_spec(F[:128],np.linspace(0,T[delta],delta))
            # flip the freq axis so lower freq is at bottom
            spec_intp = np.flip(spec_intp, axis=0).astype(np.float32)
            spec_intp = torch.from_numpy(spec_intp.copy())
            spec_wins.append(spec_intp)
            # also record the id of the window
            id_this = f'{cond_name}#{os.path.basename(fn_wav)}_{label[n]}_{n:05}_{t:05}'
            id_wins.append(id_this)
    # stack to 3d tensor
    if len(spec_wins)>0:
        spec_wins = torch.stack(spec_wins)
        return(spec_wins, id_wins)
    

def ZZ_chopSpectWinsSave_v2(fn_pair):
    # chop spectrograms into windows
    # fn_pair=[fn_wav, param_chop]
    # assume that XXX.time.txt and XXX.label.txt exist for fn_wav in the same folder
    # return a list of spectrogram windows and associated ids for each window
    [fn_wav, param_chop] = fn_pair
    cond_name = fn_wav.split('/')[-2]
    delta = param_chop['delta']  # this is the width of spectrogram window, in unit of spectrogram frames, Han uses 32
    hop = param_chop['hop']  # default hop length is 1 frame, convert to sec: hop * dt, Han uses 1
    param_spec = {'nperseg': param_chop['nperseg'], 'noverlap':param_chop['noverlap']}
    # load wav file
    fs, data = wavfile.read(fn_wav)
    data = data / 32767  # convert to float
    # band pass raw data
    data = butter_bandpass_filter(data, param_chop['bp_low'], param_chop['bp_high'], fs, order=2)
    # get spectrograms using Han's method
    [spec_full, dt, F ,T] = Han_VAE_scriptNoPitch.get_spec(data, fs, param_spec)
    # load the segmentation
    seg = np.loadtxt(fn_wav.replace('wav', 'time.txt'), delimiter=',').astype('int')  # load the segmentation
    seg = np.round(seg/fs/dt).astype('int')  # convert syllable boundary to spectrogram frames
    seg[seg>=len(T)]= len(T)
    # load the syllable labels
    label = np.loadtxt(fn_wav.replace('wav', 'label.txt'), dtype='str',delimiter=',')
    # when only have one syllable
    if len(seg.shape)==1:
        seg = np.expand_dims(seg, axis=0)
        label= np.expand_dims(label, axis=0)
    # loop through syllables, chop into spectrogram windows
    spec_wins = []
    id_wins = []
    # save the spectrograms in the same folder
    fd_save_temp = fn_wav.replace('.wav', '')
    if os.path.exists(fd_save_temp):
        shutil.rmtree(fd_save_temp)
    os.makedirs(fd_save_temp)
    for n in range(seg.shape[0]):
        spec_tmp = spec_full[:,seg[n,0]:seg[n,1]]
        # save the spectrogram of the syllables
        fn_h5 = os.path.join(fd_save_temp, f'{cond_name}#{os.path.basename(fn_wav)}_{label[n]}_{n:05}.h5')
        # print(fn_h5)
        with h5py.File(fn_h5, 'w') as hf:
            hf.create_dataset('spec', data=spec_tmp)
            hf.create_dataset('seg', data=seg[n,:])
            hf.create_dataset('seg_sec', data=seg[n,:]*dt)
        for t in range(0, spec_tmp.shape[1]-delta, hop):
            interp_spec = scipy.interpolate.RectBivariateSpline(F[:128],np.linspace(0,T[delta],delta), spec_tmp[:128,t:t+delta])
            spec_intp = interp_spec(F[:128],np.linspace(0,T[delta],delta))
            # flip the freq axis so lower freq is at bottom
            spec_intp = np.flip(spec_intp, axis=0).astype(np.float32)
            spec_intp = torch.from_numpy(spec_intp.copy())
            spec_wins.append(spec_intp)
            # also record the id of the window
            id_this = f'{cond_name}#{os.path.basename(fn_wav)}_{label[n]}_{n:05}_{t:05}'
            id_wins.append(id_this)
    # stack to 3d tensor
    if len(spec_wins)>0:
        spec_wins = torch.stack(spec_wins)
        return(spec_wins, id_wins)

def ZZ_chopSpectWinsSave_v4(fn_pair):
    # chop spectrograms into windows
    # fn_pair=[fn_wav, param_chop, fd_save_win_this]
    # assume that XXX.time.txt and XXX.label.txt exist for fn_wav in the same folder
    # return a list of spectrogram windows and associated ids for each window
    [fn_wav, param_chop, fd_save_win_this] = fn_pair
    cond_name = fn_wav.split('/')[-2]
    delta = param_chop['delta']  # this is the width of spectrogram window, in unit of spectrogram frames, Han uses 32
    hop = param_chop['hop']  # default hop length is 1 frame, convert to sec: hop * dt, Han uses 1
    param_spec = {'nperseg': param_chop['nperseg'], 'noverlap':param_chop['noverlap']}
    # load wav file
    fs, data = wavfile.read(fn_wav)
    data = data / 32767  # convert to float
    # band pass raw data
    data = butter_bandpass_filter(data, param_chop['bp_low'], param_chop['bp_high'], fs, order=2)
    # get spectrograms using Han's method
    [spec_full, dt, F ,T] = Han_VAE_scriptNoPitch.get_spec(data, fs, param_spec)
    # load the segmentation
    seg = np.loadtxt(fn_wav.replace('wav', 'time.txt'), delimiter=',').astype('int')  # load the segmentation
    seg = np.round(seg/fs/dt).astype('int')  # convert syllable boundary to spectrogram frames
    seg[seg>=len(T)]= len(T)
    # load the syllable labels
    label = np.loadtxt(fn_wav.replace('wav', 'label.txt'), dtype='str',delimiter=',')
    # when only have one syllable
    if len(seg.shape)==1:
        seg = np.expand_dims(seg, axis=0)
        label= np.expand_dims(label, axis=0)
    # loop through syllables, chop into spectrogram windows
    spec_wins = []
    id_wins = []
    # save the spectrograms in the same folder
    fd_save_temp = os.path.join(fd_save_win_this, os.path.basename(fn_wav).replace('.wav', ''))
    if os.path.exists(fd_save_temp):
        shutil.rmtree(fd_save_temp)
    os.makedirs(fd_save_temp)
    for n in range(seg.shape[0]):
        spec_tmp = spec_full[:,seg[n,0]:seg[n,1]]
        # save the spectrogram of the syllables
        fn_h5 = os.path.join(fd_save_temp, f'{cond_name}#{os.path.basename(fn_wav)}_{label[n]}_{n:05}.h5')
        # print(fn_h5)
        with h5py.File(fn_h5, 'w') as hf:
            hf.create_dataset('spec', data=spec_tmp)
            hf.create_dataset('seg', data=seg[n,:])
            hf.create_dataset('seg_sec', data=seg[n,:]*dt)
        for t in range(0, spec_tmp.shape[1]-delta, hop):
            interp_spec = scipy.interpolate.RectBivariateSpline(F[:128],np.linspace(0,T[delta],delta), spec_tmp[:128,t:t+delta])
            spec_intp = interp_spec(F[:128],np.linspace(0,T[delta],delta))
            # flip the freq axis so lower freq is at bottom
            spec_intp = np.flip(spec_intp, axis=0).astype(np.float32)
            spec_intp = torch.from_numpy(spec_intp.copy())
            spec_wins.append(spec_intp)
            # also record the id of the window
            id_this = f'{cond_name}#{os.path.basename(fn_wav)}_{label[n]}_{n:05}_{t:05}'
            id_wins.append(id_this)
    # stack to 3d tensor
    if len(spec_wins)>0:
        spec_wins = torch.stack(spec_wins)
        return(spec_wins, id_wins)
    
def ZZ_chopSpectWinsSave_v5(fn_pair):
    # chop spectrograms into windows
    # fn_pair=[fn_wav, param_chop, fd_save_win_this]
    # assume that XXX.time.txt and XXX.label.txt exist for fn_wav in the same folder
    # return a list of spectrogram windows and associated ids for each window
    # modified: ignore ultra short segments
    [fn_wav, param_chop, fd_save_win_this] = fn_pair
    cond_name = fn_wav.split('/')[-2]
    delta = param_chop['delta']  # this is the width of spectrogram window, in unit of spectrogram frames, Han uses 32
    hop = param_chop['hop']  # default hop length is 1 frame, convert to sec: hop * dt, Han uses 1
    param_spec = {'nperseg': param_chop['nperseg'], 'noverlap':param_chop['noverlap']}
    # load wav file
    fs, data = wavfile.read(fn_wav)
    if len(data.shape)>1:
        data = data[:, param_chop['focal_ch']]
    data = data / 32767  # convert to float
    # band pass raw data
    data = butter_bandpass_filter(data, param_chop['bp_low'], param_chop['bp_high'], fs, order=2)
    # get spectrograms using Han's method
    [spec_full, dt, F ,T] = Han_VAE_scriptNoPitch.get_spec(data, fs, param_spec)
    # load the segmentation
    seg = np.loadtxt(fn_wav.replace('wav', 'time.txt'), delimiter=',').astype('int')  # load the segmentation
    seg = np.round(seg/fs/dt).astype('int')  # convert syllable boundary to spectrogram frames
    seg[seg>=len(T)]= len(T)
    # load the syllable labels
    label = np.loadtxt(fn_wav.replace('wav', 'label.txt'), dtype='str',delimiter=',')
    # when only have one syllable
    if len(seg.shape)==1:
        seg = np.expand_dims(seg, axis=0)
        label= np.expand_dims(label, axis=0)
    # loop through syllables, chop into spectrogram windows
    spec_wins = []
    id_wins = []
    # save the spectrograms in the same folder
    fd_save_temp = os.path.join(fd_save_win_this, os.path.basename(fn_wav).replace('.wav', ''))
    if os.path.exists(fd_save_temp):
        shutil.rmtree(fd_save_temp)
    os.makedirs(fd_save_temp)
    for n in range(seg.shape[0]):
        spec_tmp = spec_full[:,seg[n,0]:seg[n,1]]
        # save the spectrogram of the syllables
        fn_h5 = os.path.join(fd_save_temp, f'{cond_name}#{os.path.basename(fn_wav)}_{label[n]}_{n:05}.h5')
        # print(fn_h5)
        with h5py.File(fn_h5, 'w') as hf:
            hf.create_dataset('spec', data=spec_tmp)
            hf.create_dataset('seg', data=seg[n,:])
            hf.create_dataset('seg_sec', data=seg[n,:]*dt)
        for t in range(0, spec_tmp.shape[1]-delta, hop):
            interp_spec = scipy.interpolate.RectBivariateSpline(F[:128],np.linspace(0,T[delta],delta), spec_tmp[:128,t:t+delta])
            spec_intp = interp_spec(F[:128],np.linspace(0,T[delta],delta))
            # flip the freq axis so lower freq is at bottom
            spec_intp = np.flip(spec_intp, axis=0).astype(np.float32)
            spec_intp = torch.from_numpy(spec_intp.copy())
            spec_wins.append(spec_intp)
            # also record the id of the window
            id_this = f'{cond_name}#{os.path.basename(fn_wav)}_{label[n]}_{n:05}_{t:05}'
            id_wins.append(id_this)
    # stack to 3d tensor
    if len(spec_wins)>0:
        spec_wins = torch.stack(spec_wins)
        return(spec_wins, id_wins)
    
def ZZ_chopSpectWinsSave_v6(fn_pair):
    # chop spectrograms into windows
    # fn_pair=[fn_wav, param_chop, fd_save_win_this]
    # assume that XXX.time.txt and XXX.label.txt exist for fn_wav in the same folder
    # return a list of spectrogram windows and associated ids for each window
    # modified: ignore ultra short segments
    [fn_wav, param_chop, fd_save_win_this] = fn_pair
    delta = param_chop['delta']  # this is the width of spectrogram window, in unit of spectrogram frames, Han uses 32
    hop = param_chop['hop']  # default hop length is 1 frame, convert to sec: hop * dt, Han uses 1
    param_spec = {'nperseg': param_chop['nperseg'], 'noverlap':param_chop['noverlap']}
    # load wav file
    fs, data = wavfile.read(fn_wav)
    if len(data.shape)>1:
        data = data[:, param_chop['focal_ch']]
    data = data / 32767  # convert to float
    # band pass raw data
    data = butter_bandpass_filter(data, param_chop['bp_low'], param_chop['bp_high'], fs, order=2)
    # get spectrograms using Han's method
    [spec_full, dt, F ,T] = Han_VAE_scriptNoPitch.get_spec(data, fs, param_spec)
    # load the segmentation
    seg = np.loadtxt(fn_wav.replace('wav', 'time.txt'), delimiter=',').astype('int')  # load the segmentation
    seg = np.round(seg/fs/dt).astype('int')  # convert syllable boundary to spectrogram frames
    seg[seg>=len(T)]= len(T)
    # load the syllable labels
    label = np.loadtxt(fn_wav.replace('wav', 'label.txt'), dtype='str',delimiter=',')
    # when only have one syllable
    if len(seg.shape)==1:
        seg = np.expand_dims(seg, axis=0)
        label= np.expand_dims(label, axis=0)
    # loop through syllables, chop into spectrogram windows
    spec_wins = []
    id_wins = []
    fns_wav_win = []  # keep a note of the original wav file
    # save the spectrograms in the same folder
    fd_save_temp = os.path.join(fd_save_win_this, os.path.basename(fn_wav).replace('.wav', ''))
    if os.path.exists(fd_save_temp):
        shutil.rmtree(fd_save_temp)
    os.makedirs(fd_save_temp)
    for n in range(seg.shape[0]):
        spec_tmp = spec_full[:,seg[n,0]:seg[n,1]]
        # save the spectrogram of the syllables
        fn_h5 = os.path.join(fd_save_temp, f'{os.path.basename(fn_wav)}_{label[n]}_{n:05}.h5')
        # print(fn_h5)
        with h5py.File(fn_h5, 'w') as hf:
            hf.create_dataset('spec', data=spec_tmp)
            hf.create_dataset('seg', data=seg[n,:])
            hf.create_dataset('seg_sec', data=seg[n,:]*dt)
        for t in range(0, spec_tmp.shape[1]-delta, hop):
            interp_spec = scipy.interpolate.RectBivariateSpline(F[:128],np.linspace(0,T[delta],delta), spec_tmp[:128,t:t+delta])
            spec_intp = interp_spec(F[:128],np.linspace(0,T[delta],delta))
            # flip the freq axis so lower freq is at bottom
            spec_intp = np.flip(spec_intp, axis=0).astype(np.float32)
            spec_intp = torch.from_numpy(spec_intp.copy())
            spec_wins.append(spec_intp)
            # also record the id of the window
            id_this = f'{os.path.basename(fn_wav)}_{label[n]}_{n:05}_{t:05}'
            id_wins.append(id_this)
            fns_wav_win.append(fn_wav)
    # stack to 3d tensor
    if len(spec_wins)>0:
        spec_wins = torch.stack(spec_wins)
        return(spec_wins, id_wins, fns_wav_win)
    
def ZZ_chopSpectWinsSave_v7(fn_pair):
    # chop spectrograms into windows
    # fn_pair=[fn_wav, param_chop, fd_save_win_this]
    # assume that XXX.time.txt and XXX.label.txt exist for fn_wav in the same folder
    # return a list of spectrogram windows and associated ids for each window
    # modified: ignore ultra short segments
    [fn_wav, param_chop, fd_save_win_this] = fn_pair
    delta = param_chop['delta']  # this is the width of spectrogram window, in unit of spectrogram frames, Han uses 32
    hop = param_chop['hop']  # default hop length is 1 frame, convert to sec: hop * dt, Han uses 1
    param_spec = {'nperseg': param_chop['nperseg'], 'noverlap':param_chop['noverlap']}
    # load wav file
    fs, data = wavfile.read(fn_wav)
    if len(data.shape)>1:
        data = data[:, param_chop['focal_ch']]
    data = data / 32767  # convert to float
    # band pass raw data
    data = butter_bandpass_filter(data, param_chop['bp_low'], param_chop['bp_high'], fs, order=2)
    # get spectrograms using Han's method
    [spec_full, dt, F ,T] = Han_VAE_scriptNoPitch.get_spec(data, fs, param_spec)
    # load the segmentation
    seg = np.loadtxt(fn_wav.replace('wav', 'time.txt'), delimiter=',').astype('int')  # load the segmentation
    seg = np.round(seg/fs/dt).astype('int')  # convert syllable boundary to spectrogram frames
    seg[seg>=len(T)]= len(T)
    # load the syllable labels
    label = np.loadtxt(fn_wav.replace('wav', 'label.txt'), dtype='str',delimiter=',', encoding='utf-8-sig')
    # when only have one syllable
    if len(seg.shape)==1:
        seg = np.expand_dims(seg, axis=0)
        label= np.expand_dims(label, axis=0)
    # loop through syllables, chop into spectrogram windows
    spec_wins = []
    id_wins = []
    fns_wav_win = []  # keep a note of the original wav file
    # save the spectrograms in the same folder
    fd_save_temp = os.path.join(fd_save_win_this, os.path.basename(fn_wav).replace('.wav', ''))
    if os.path.exists(fd_save_temp):
        shutil.rmtree(fd_save_temp)
    os.makedirs(fd_save_temp)
    for n in range(seg.shape[0]):
        spec_tmp = spec_full[:,seg[n,0]:seg[n,1]]
        # save the spectrogram of the syllables
        fn_h5 = os.path.join(fd_save_temp, f'{os.path.basename(fn_wav)}_{label[n]}_{n:05}.h5')
        # print(fn_h5)
        with h5py.File(fn_h5, 'w') as hf:
            hf.create_dataset('spec', data=spec_tmp)
            hf.create_dataset('seg', data=seg[n,:])
            hf.create_dataset('seg_sec', data=seg[n,:]*dt)
        for t in range(0, spec_tmp.shape[1]-delta, hop):
            interp_spec = scipy.interpolate.RectBivariateSpline(F[:128],np.linspace(0,T[delta],delta), spec_tmp[:128,t:t+delta])
            spec_intp = interp_spec(F[:128],np.linspace(0,T[delta],delta))
            # flip the freq axis so lower freq is at bottom
            spec_intp = np.flip(spec_intp, axis=0).astype(np.float32)
            spec_intp = torch.from_numpy(spec_intp.copy())
            spec_wins.append(spec_intp)
            # also record the id of the window
            id_this = f'{os.path.basename(fn_wav)}_{label[n]}_{n:05}_{t:05}'
            id_wins.append(id_this)
            fns_wav_win.append(fn_wav)
    # stack to 3d tensor
    if len(spec_wins)>0:
        spec_wins = torch.stack(spec_wins)
        return(spec_wins, id_wins, fns_wav_win)
    
class DatasetSpectWinsVAE(Dataset):
    # Dataset class for loading spectrogram windows during training
    def __init__(self, fns_h5):
        self.fns_h5 = fns_h5
        
    def __len__(self):
        return(len(self.fns_h5))
    
    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        # read h5 file
        with h5py.File(self.fns_h5[idx], 'r') as hf:
            spec = np.array(hf.get('spec'), dtype='float32')
        return spec
    

def vae_data_loaders(fns_h5,split=0.8,batch_size=10, shuffle=True, num_workers=4):
    # split the dataset into train and split, return a dict of dataloaders
    assert(split > 0.0 and split <= 1.0)
    if shuffle:
        np.random.seed(42)
        perm = np.random.permutation(len(fns_h5))
        fns_h5 = [fns_h5[aa] for aa in perm]
        np.random.seed(None)
    # Split
    i = int(round(split * len(fns_h5)))
    fns_train = fns_h5[:i]
    fns_test = fns_h5[i:]

    # construct the dataset and dataloader
    ds_train = DatasetSpectWinsVAE(fns_train)
    train_dataloader = DataLoader(ds_train, batch_size=batch_size, shuffle=shuffle, num_workers=num_workers)
    ds_test = DatasetSpectWinsVAE(fns_test)
    test_dataloader = DataLoader(ds_test, batch_size=batch_size, shuffle=shuffle, num_workers=num_workers)
    return {'train':train_dataloader, 'test':test_dataloader}


def ZZ_readNC(fn):
    """a function to read netcdf (.nc) files"""
    with netcdf_file(fn, 'r', mmap=False) as nf:
        dt = nf.variables['dt'].data
        fs = int(1/dt)
        d = nf.variables['data'][:].copy()
    return(d, fs)


def ZZ_makeSpectrogramSyllable_v1(fn_pair):
    # a function to generate spectrograms from labeled syllables, then embed it in a predefined matrix
    fn_sound, onsets, offsets, embed_width, fd_save_win_this, target, param = fn_pair
    # read the sound file, could be nc or wav format
    sound_format = fn_sound.split('.')[-1]
    if sound_format=='nc':
        data, fs = ZZ_readNC(fn_sound)
    elif sound_format=='wav':
        fs, data = wavfile.read(fn_sound)
    else:
        print("Sound file format needs to be wav or nc (netcdf)!")
    # bandpass
    data = butter_bandpass_filter(data, param['bp_low'], param['bp_high'], fs, order=2)
    # get spectrograms using Han's method
    [spec_full, dt, F ,T] = Han_VAE_scriptNoPitch.get_spec(data, fs, param)
    # when the power at certain point is below certain threshold, set it to that threshold
    spec_full[np.where(spec_full<param['mask_thre'])] = param['mask_thre']
    # loop through the annotated syllables
    spec_all = np.zeros([len(onsets), spec_full.shape[0], int(embed_width/dt)]) + param['mask_thre']
    for ii in range(len(onsets)):
        # convert onset/offset time (sec) to frame index in spectrograms
        idx_onset = int(onsets[ii]/dt)
        idx_offset = int(offsets[ii]/dt)
        spec_this = spec_full[:, idx_onset:idx_offset]
        # flip the y-axis so low-frequency is at the bottom
        spec_this = np.flip(spec_this, axis=0)
        # calculate where to place the spectrograms in the embed matrix
        pos_off = int((embed_width-(offsets[ii]-onsets[ii]))/dt/2)
        if pos_off>=0:  # syllable is shorter than embed matrix
            spec_all[ii, :, pos_off:(pos_off+spec_this.shape[1])] = spec_this
        else:  # syllable is longer than embed matrix, keep the center part
            spec_all[ii,:,0:spec_all.shape[2]] = spec_this[:,(-pos_off):(-pos_off+spec_all.shape[2])]
    # also save the spectrograms as a h5 file for future visualization
    fd_save_temp = os.path.join(fd_save_win_this, os.path.basename(fn_sound).replace(f'.{sound_format}', ''))
    if os.path.exists(fd_save_temp):
        shutil.rmtree(fd_save_temp)
    os.makedirs(fd_save_temp)
    cond_name = fd_save_temp.split('/')[-2]
    ids_wins = []
    for ii in range(spec_all.shape[0]):
        fn_h5 = os.path.join(fd_save_temp, f'{cond_name}#{os.path.basename(fn_sound)}_{target}_{ii:05}.h5')
        ids_wins.append(f'{cond_name}#{os.path.basename(fn_sound)}_{target}_{ii:05}')
        with h5py.File(fn_h5, 'w') as hf:
            hf.create_dataset('spec', data=spec_all[ii,:,:])
            hf.create_dataset('seg', data=np.array([int(onsets[ii]/dt), int(offsets[ii]/dt)]))
            hf.create_dataset('seg_sec', data=np.array([onsets[ii], offsets[ii]]))
    # only take the first 128 dims on frequency
    spec_all = spec_all[:, 0:128,:]
    spec_all = torch.Tensor(spec_all.astype(np.float32))
    return(spec_all, ids_wins)


def getAudioSpectrogramZZ_flexible_v1(audio, fs, NFFT=512, windowSize=512, hopsec=0.001, flim=[500,8000], clim=[-3,7]):
    # a function to generate spectrograms using Brian's electro_gui script
    windowOverlap = windowSize-np.floor(hopsec*fs);
    # [S,F,t] = spectrogram(audio, windowSize, windowOverlap, NFFT, fs);
    F,t,S = scipy.signal.stft(audio, fs, window='hann', nperseg=windowSize, noverlap=windowOverlap, nfft=NFFT)
    dt = t[1] - t[0]
    # only keep the meaningful frequency range
    freqInRange = np.where((F>=flim[0]) & (F<=flim[1]))[0]
    f = F[freqInRange]
    # calculate the power spectrogram
    power = 2*np.log(np.abs(S[freqInRange,:])+1e-16)+20  # offset by 20db
    # save the original power before coloring
    powerRaw = power.copy(); 
    # coloring the spectrogram based on clim
    cmap_temp = plt.cm.jet; 
    c = cmap_temp(np.linspace(0, 1, 1024))
    numColors = c.shape[0]
    c[0, :] = [0, 0, 0, 1]
    # set values outside clim range to the boundary condition
    # power = np.floor((numColors - 1) * (power - clim[0]) / (clim[1] - clim[0]) + 1)
    power = (numColors - 1) * (power - clim[0]) / (clim[1] - clim[0]) + 1
    power[np.where(power < 0)] = 0
    power[np.where(power > numColors-1)] = numColors-1
    powerGrey = power.copy()
    power = np.floor(power)
    power = power.astype(np.int64)
    power = c[power]
    return (power, powerRaw, powerGrey, S, f, t, dt)


def ZZ_interp_spec(original_array, desired_dim=[128, 512]):
    # a function to interpolate spectrogram into arbitrary size
    x = np.linspace(0, 1, original_array.shape[0])
    y = np.linspace(0, 1, original_array.shape[1])
    # New indices for the target size 128x512
    x_new = np.linspace(0, 1, desired_dim[0])
    y_new = np.linspace(0, 1, desired_dim[1])
    # Create a bivariate spline interpolator
    interpolator = scipy.interpolate.RectBivariateSpline(x, y, original_array)
    # Evaluate the spline on the new grid
    resized_array = interpolator(x_new, y_new)
    return resized_array


def ZZ_chopSpectWinsSave_v3(fn_pair):
    # chop spectrograms into windows
    # fn_pair=[fn_wav, param_chop]
    # assume that XXX.time.txt and XXX.label.txt exist for fn_wav in the same folder
    # return a list of spectrogram windows and associated ids for each window
    [fn_wav, param_chop] = fn_pair
    cond_name = fn_wav.split('/')[-2]
    delta = param_chop['chopWinColumns']  # this is the width of spectrogram window, in unit of spectrogram frames, Han uses 32
    hop = param_chop['chopHopColumns']  # default hop length is 1 frame, convert to sec: hop * dt, Han uses 1
    # load wav file
    # fs, data = wavfile.read(fn_wav)
    # data = data / 32767  # convert to float
    data, fs = librosa.load(fn_wav, sr=None)
    # band pass raw data
    data = butter_bandpass_filter(data, param_chop['bp_low'], param_chop['bp_high'], fs, order=2)
    # get spectrograms using Han's method
    power, powerRaw, powerGrey, S, f, T, dt = getAudioSpectrogramZZ_flexible_v1(data, fs, NFFT=param_chop['NFFT'], windowSize=param_chop['specWinSize'], 
                                                                            hopsec=param_chop['specHopSec'], flim=param_chop['flim'], clim=param_chop['clim'])
    # load the segmentation
    seg = np.loadtxt(fn_wav.replace('wav', 'time.txt'), delimiter=',').astype('int')  # load the segmentation
    seg = np.round(seg/fs/dt).astype('int')  # convert syllable boundary to spectrogram frames
    seg[seg>=len(T)]= len(T)
    #load the syllable labels
    label = np.loadtxt(fn_wav.replace('wav', 'label.txt'), dtype='str',delimiter=',')
    # when only have one syllable
    if len(seg.shape)==1:
        seg = np.expand_dims(seg, axis=0)
        label= np.expand_dims(label, axis=0)
    # loop through syllables, chop into spectrogram windows
    spec_wins = []
    id_wins = []
    # save the spectrograms in the same folder
    fd_save_temp = fn_wav.replace('.wav', '')
    if os.path.exists(fd_save_temp):
        shutil.rmtree(fd_save_temp)
    os.makedirs(fd_save_temp)
    for n in range(seg.shape[0]):
        # n = 1
        spec_tmp = powerGrey[:,seg[n,0]:seg[n,1]]
        # save the spectrogram of the syllables
        fn_h5 = os.path.join(fd_save_temp, f'{cond_name}#{os.path.basename(fn_wav)}_{label[n]}_{n:05}.h5')
        # print(fn_h5)
        with h5py.File(fn_h5, 'w') as hf:
            hf.create_dataset('spec', data=spec_tmp)
            hf.create_dataset('seg', data=seg[n,:])
            hf.create_dataset('seg_sec', data=seg[n,:]*dt)
            for t in range(0, spec_tmp.shape[1]-delta, hop):
                spec_this = spec_tmp[:, t:(t+delta)]
                # reshape the spectrogram window to define size
                spec_intp = ZZ_interp_spec(spec_this, [128, 512])
                # remove arbitrary noise introduced in intepolation
                spec_intp[spec_intp<0] = 0
                # flip the freq axis so lower freq is at bottom
                spec_intp = np.flip(spec_intp, axis=0).astype(np.float32)
                spec_intp = torch.from_numpy(spec_intp.copy())
                spec_wins.append(spec_intp)
                # also record the id of the window
                id_this = f'{cond_name}#{os.path.basename(fn_wav)}_{label[n]}_{n:05}_{t:05}'
                id_wins.append(id_this)
    # stack to 3d tensor
    if len(spec_wins)>0:
        spec_wins = torch.stack(spec_wins)
        return(spec_wins, id_wins)
    
    
def SpecSyllableSave_v3(fn, params, thres_t):
    # given a wav file and its segment time/label files, calculate the spectrograms
    # save syllable spectrograms into a separate folder
    # load the wav signal and associated segmentations
    wav_tmp, fs = librosa.load(fn, sr=None)
    segments = np.loadtxt(fn.replace('.wav', '.time.txt'), delimiter=',', dtype='int')
    labels = np.loadtxt(fn.replace('.wav', '.label.txt'), delimiter=',', dtype='str')
    # ignore syllables that are shorter than thres
    thres_p = int(thres_t*fs)  # convert to the number of data points
    # dump the spectrograms of syllables into pickle files in a separate folder
    fd_save_this = fn.replace('.wav', '')
    if os.path.exists(fd_save_this):
        shutil.rmtree(fd_save_this)
    os.makedirs(fd_save_this)
    # define a within-scope function to save spectrograms of syllables
    def saveSegSpec(m):
         # calculate spectrograms of the segment only
        wav_this = wav_tmp[segments[m,0]:segments[m,1]]
        if wav_this.shape[0]>=thres_p:
            power, powerRaw, powerGrey, S, f, t, dt = getAudioSpectrogramZZ_flexible_v1(audio=wav_this, fs=fs, NFFT=params['NFFT'], 
                                                    windowSize=params['windowSize'], hopsec=params['hopsec'], flim=params['flim'], clim=params['clim'])
            tmp_data = {}
            file_save_name = f'{os.path.basename(fn)}_seg{m:08d}.p'
            tmp_data['filename'] = fn
            # tmp_data['budgieID'] = budgieID
            tmp_data['syl_id'] = m
            tmp_data['label'] = labels[m]
            tmp_data['spec'] = powerGrey/1024.0  # convert values to between 0 and 1
            tmp_data['dt'] = dt
            tmp_data['F'] = f
            tmp_data['start_idx'] = segments[m,0]
            tmp_data['end_idx'] = segments[m,1]
            pickle.dump(tmp_data, open(os.path.join(fd_save_this, file_save_name), 'wb'))
    with Parallel(n_jobs=num_cpu, verbose=0) as parallel:
        parallel(delayed(saveSegSpec)(m) for m in range(segments.shape[0]))
    return fd_save_this


def intp_spec_v2(original_array, target_shape=[128,32]):
    # interpret spectrogram into a target size
    # Define the dimensions of the original and target arrays
    ori_F = np.linspace(0, original_array.shape[0]-1, original_array.shape[0])
    ori_T = np.linspace(0, original_array.shape[1]-1, original_array.shape[1])
    interp=RegularGridInterpolator([ori_T, ori_F], original_array.T)
    target_F = np.linspace(0, original_array.shape[0]-1, target_shape[0])
    target_T = np.linspace(0, original_array.shape[1]-1, target_shape[1])
    xx =np.meshgrid(target_T, target_F, indexing='ij')
    sample_points = np.array([xx[0].ravel(), xx[1].ravel()]).T 
    resized_array = interp(sample_points, method='linear').reshape(target_shape[1],target_shape[0]).T
    return(resized_array)
    
    
def ZZ_createDataset_v2(spec_dirs, fn_h5, chop_frames, sliding_frames, spec_input):
    # a function to create training dataset from a list of syllable spectrogram files
    # first count the number of total windows so h5 can be initated
    def countChopWin_v2(fn):
        # distinguish two cases: 
        # 1. for syllables longer than chop size, count how many spectrogram chops
        # 2. for syllables shorter than chop size, will embed it in the middle, count 1 
        with open(fn, 'rb') as file:
            tmp_data = pickle.load(file)
            spec = tmp_data['spec']
        # how many spectrogram chops?
        if spec.shape[1]>chop_frames:
            chop_start = np.arange(0, spec.shape[1]-chop_frames, sliding_frames)
            return(chop_start.shape[0])
        else:
            return(1)
    with Parallel(n_jobs=num_cpu, verbose=1) as parallel:
        results = parallel(delayed(countChopWin_v2)(fn) for fn in spec_dirs)
    img_cnt = np.sum(results)
    print(f'total number of spec windows: {img_cnt}')
    # then create the h5 file
    with h5py.File(fn_h5, 'w') as hf:
        hf.create_dataset('spec_data', shape=(img_cnt, spec_input[0], spec_input[1]))
        count = 0
        id_pd = pd.DataFrame()  # save the id of each spectrogram chop in a dataframe
        for n in tqdm.tqdm(range(len(spec_dirs))):
        # for n in tqdm.tqdm(range(2)):
            fn = spec_dirs[n]
            with open(fn, 'rb') as file:
                tmp_data = pickle.load(file)
                spec = tmp_data['spec']
            if spec.shape[1]>chop_frames:  # syllable is longer than chop size: chop
                chop_start = np.arange(0, spec.shape[1]-chop_frames, sliding_frames)
                for i in range(len(chop_start)):
                    spec_this = spec[:, chop_start[i]:(chop_start[i]+chop_frames)]
                    # resize the spec
                    spec_this = intp_spec_v2(spec_this, target_shape=spec_input)
                    hf['spec_data'][count,:,:] = spec_this
                    id_this = pd.DataFrame([{'fn_wav': tmp_data['filename'], 'label':tmp_data['label'], 'fn_p':fn, 'idx_start':chop_start[i], 'idx_end': chop_start[i]+chop_frames}])
                    id_pd = pd.concat([id_pd, id_this], ignore_index=True)
                    count += 1
            else: # syllable is shorter than chop size: embed in the middle
                temp = np.zeros((spec.shape[0], chop_frames))
                off_idx = int((chop_frames-spec.shape[1])/2)
                temp[:,off_idx:(off_idx+spec.shape[1])] = spec
                # resize the spec
                spec_this = intp_spec_v2(temp, target_shape=spec_input)
                hf['spec_data'][count,:,:] = spec_this
                id_this = pd.DataFrame([{'fn_wav': tmp_data['filename'], 'label':tmp_data['label'], 'fn_p':fn, 'idx_start':0, 'idx_end': spec.shape[1]}])
                id_pd = pd.concat([id_pd, id_this], ignore_index=True)
                count += 1

    return(count, id_pd)


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


def ZZ_chopSpectWinsSave_v8(fn_pair):
    # chop spectrograms into windows
    # fn_pair=[fn_wav, param_chop, fd_save_win_this]
    # assume that XXX.time.txt and XXX.label.txt exist for fn_wav in the same folder
    # return a list of spectrogram windows and associated ids for each window
    [fn_wav, param_chop, fd_save_win_this] = fn_pair
    delta = param_chop['delta']  # this is the width of spectrogram window, in unit of spectrogram frames, Han uses 32
    hop = param_chop['hop']  # default hop length is 1 frame, convert to sec: hop * dt, Han uses 1
    # load wav file
    fs, data = wavfile.read(fn_wav)
    if len(data.shape)>1:
        data = data[:, param_chop['focal_ch']]
    data = data / 32767  # convert to float
    # band pass raw audio data
    data = butter_bandpass_filter(data, param_chop['target_freqs'][0], param_chop['target_freqs'][-1], fs, order=2)
    # load the segmentation
    seg = np.loadtxt(fn_wav.replace('wav', 'time.txt'), delimiter=',').astype('int')  # load the segmentation
    # load the syllable labels
    label = np.loadtxt(fn_wav.replace('wav', 'label.txt'), dtype='str',delimiter=',', encoding='utf-8-sig')
    # when only have one syllable
    if len(seg.shape)==1:
        seg = np.expand_dims(seg, axis=0)
        label= np.expand_dims(label, axis=0)
    # loop through each segment, calculate spectrogram, then chop
    spec_wins = []
    id_wins = []
    fns_wav_win = []  # keep a note of the original wav file
    # save the spectrograms in the same folder
    fd_save_temp = os.path.join(fd_save_win_this, os.path.basename(fn_wav).replace('.wav', ''))
    if os.path.exists(fd_save_temp):
        shutil.rmtree(fd_save_temp)
    os.makedirs(fd_save_temp)
    for n in range(seg.shape[0]):
        if (seg[n,1]-seg[n,0]) < param_chop['nperseg']:  # ignore ultral short segment
            continue              
        spec_tmp, dt, f, t = ZZ_get_spec(seg[n,0], seg[n,1], data, fs, nperseg=param_chop['nperseg'], noverlap=param_chop['noverlap'], 
                                                   target_freqs=param_chop['target_freqs'], target_time_dim=param_chop['target_time_dim'],
                                                   bool_clip=param_chop['bool_clip'], clim=param_chop['clim'])
        # save the spectrogram of the entire syllable
        fn_h5 = os.path.join(fd_save_temp, f'{os.path.basename(fn_wav)}_{label[n]}_{n:06}.h5')
        # print(fn_h5)
        with h5py.File(fn_h5, 'w') as hf:
            hf.create_dataset('spec', data=spec_tmp)
            hf.create_dataset('seg', data=seg[n,:])
            hf.create_dataset('seg_sec', data=seg[n,:]*dt)
        # chop syllable into spectrogram windows, load into memory
        for t in range(0, spec_tmp.shape[1]-delta, hop):
            spec_intp = spec_tmp[:,t:t+delta]
            # flip the freq axis so lower freq is at bottom
            spec_intp = np.flip(spec_intp, axis=0).astype(np.float32)
            spec_intp = torch.from_numpy(spec_intp.copy())
            spec_wins.append(spec_intp)
            # also record the id of the window
            id_this = f'{os.path.basename(fn_wav)}_{label[n]}_{n:06}_{t:05}'
            id_wins.append(id_this)
            fns_wav_win.append(fn_wav)
    # stack to 3d tensor
    if len(spec_wins)>0:
        spec_wins = torch.stack(spec_wins)
    return(spec_wins, id_wins, fns_wav_win)



def SpecSyllableSave_v4(fn, fd_save_this, param_chop):
    # given a wav file and its segment time/label files in the same wav folder, calculate the spectrograms
    # save syllable spectrograms into a separate folder
    # load the wav signal and associated segmentations
    data, fs = librosa.load(fn, sr=None)
    segments = np.loadtxt(fn.replace('.wav', '.time.txt'), delimiter=',', dtype='int')
    labels = np.loadtxt(fn.replace('.wav', '.label.txt'), delimiter=',', dtype='str', encoding='utf-8-sig')
    # ignore syllables that are shorter than one spectrogram window
    thres_p = int(param_chop['nperseg']+(param_chop['delta']-1)*(param_chop['nperseg']-param_chop['noverlap']))
    # dump the spectrograms of syllables into pickle files in a separate folder
    # define a within-scope function to save spectrograms of syllables
    def saveSegSpec(m):
         # calculate spectrograms of the segment only
        if (segments[m,1]-segments[m,0])>=thres_p:
            spec_tmp, dt, f, t = ZZ_get_spec(segments[m,0], segments[m,1], data, fs, nperseg=param_chop['nperseg'], noverlap=param_chop['noverlap'], 
                                                       target_freqs=param_chop['target_freqs'], target_time_dim=param_chop['target_time_dim'],
                                                       bool_clip=param_chop['bool_clip'], clim=param_chop['clim'])
            tmp_data = {}
            file_save_name = f'{os.path.basename(fn)}_seg{m:08d}.p'
            tmp_data['filename'] = fn
            # tmp_data['budgieID'] = budgieID
            tmp_data['syl_id'] = m
            tmp_data['label'] = labels[m]
            tmp_data['spec'] = spec_tmp 
            tmp_data['dt'] = dt
            tmp_data['F'] = f
            tmp_data['start_idx'] = segments[m,0]
            tmp_data['end_idx'] = segments[m,1]
            pickle.dump(tmp_data, open(os.path.join(fd_save_this, file_save_name), 'wb'))
    with Parallel(n_jobs=num_cpu, verbose=0) as parallel:
        parallel(delayed(saveSegSpec)(m) for m in range(segments.shape[0]))
    return fd_save_this


def ZZ_createDataset_v3(spec_dirs, fn_h5, param_chop):
# a function to create VAE dataset from a list of syllable spectrogram files
# first count the number of total windows so h5 can be initated
    def countChopWin_v2(fn):
        # distinguish two cases: 
        # 1. for syllables longer than chop size, count how many spectrogram chops
        # 2. for syllables shorter than chop size, will embed it in the middle, count 1 
        with open(fn, 'rb') as file:
            tmp_data = pickle.load(file)
            spec = tmp_data['spec']
        # how many spectrogram chops?
        if spec.shape[1]>=param_chop['delta']:
            chop_start = np.arange(0, spec.shape[1]-param_chop['delta']+1, param_chop['hop'])
            return(chop_start.shape[0])
        else:
            return(0)
    with Parallel(n_jobs=num_cpu, verbose=0) as parallel:
        results = parallel(delayed(countChopWin_v2)(fn) for fn in spec_dirs)
    img_cnt = np.sum(results)
    print(f'total number of spec windows: {img_cnt}')
    # read in one spectrogram to determine the dimension
    with open(spec_dirs[0], 'rb') as file:
        loaded_object = pickle.load(file)
        spec_input = loaded_object['spec']
    # # then create the h5 file
    with h5py.File(fn_h5, 'w') as hf:
        hf.create_dataset('spec_data', shape=(img_cnt, spec_input.shape[0], param_chop['delta']))
        count = 0
        id_pd = pd.DataFrame()  # save the id of each spectrogram chop in a dataframe
        for n in tqdm.tqdm(range(len(spec_dirs))):
            fn = spec_dirs[n]
            with open(fn, 'rb') as file:
                tmp_data = pickle.load(file)
                spec = tmp_data['spec']
            if spec.shape[1]>=param_chop['delta']:  # syllable is equal or longer than chop size: chop
                chop_start = np.arange(0, spec.shape[1]-param_chop['delta']+1, param_chop['hop'])
                for i in range(len(chop_start)):
                    spec_this = spec[:, chop_start[i]:(chop_start[i]+param_chop['delta'])]
                    # add to the large h5 file
                    hf['spec_data'][count,:,:] = spec_this
                    id_this = pd.DataFrame([{'fn_wav': tmp_data['filename'], 'label':tmp_data['label'], 'fn_p':fn, 'idx_start':chop_start[i], 'idx_end': chop_start[i]+param_chop['delta']}])
                    id_pd = pd.concat([id_pd, id_this], ignore_index=True)
                    count += 1

    return(count, id_pd)


def SpecSyllableSave_v5(fn, fd_save_this, param_chop):
    # given a wav file and its segment time/label files in the same wav folder, calculate the spectrograms
    # save syllable spectrograms into a separate folder
    # load the wav signal and associated segmentations
    data, fs = librosa.load(fn, sr=None)
    segments = np.loadtxt(fn.replace('.wav', '.time.txt'), delimiter=',', dtype='int')
    labels = np.loadtxt(fn.replace('.wav', '.label.txt'), delimiter=',', dtype='str', encoding='utf-8-sig')
    # ignore syllables that are shorter than one spectrogram window
    thres_p = int(param_chop['nperseg']+(param_chop['delta']-1)*(param_chop['nperseg']-param_chop['noverlap']))
    # dump the spectrograms of syllables into pickle files in a separate folder
    # define a within-scope function to save spectrograms of syllables
    def saveSegSpec(m):
         # calculate spectrograms of the segment only
        if (segments[m,1]-segments[m,0])>=thres_p:
            spec_tmp, dt, f, t = ZZ_get_spec(segments[m,0], segments[m,1], data, fs, nperseg=param_chop['nperseg'], noverlap=param_chop['noverlap'], 
                                                       target_freqs=param_chop['target_freqs'], target_time_dim=param_chop['target_time_dim'],
                                                       bool_clip=param_chop['bool_clip'], clim=param_chop['clim'])
            tmp_data = {}
            # save spectrogram windows as h5 file
            fn_h5 = os.path.join(fd_save_this, f'{os.path.basename(fn)}_{m:08d}_{labels[m]}.h5')
            # print(fn_h5)
            with h5py.File(fn_h5, 'w') as hf:
                hf.create_dataset('spec', data=spec_tmp)
                hf.create_dataset('fn_wav', data=fn)
                hf.create_dataset('start_idx', data=segments[m,0])
                hf.create_dataset('end_idx', data=segments[m,1])
                hf.create_dataset('fs', data=fs)
                hf.create_dataset('f', data=f)
                hf.create_dataset('t', data=t)
                hf.create_dataset('dt', data=dt)
                hf.create_dataset('filename', data=fn)
                hf.create_dataset('label', data=list(labels[m]))
                hf.create_dataset('syl_id', data=m)
    with Parallel(n_jobs=num_cpu, verbose=0) as parallel:
        parallel(delayed(saveSegSpec)(m) for m in range(segments.shape[0]))
    return fd_save_this


def ZZ_createDataset_v4(spec_dirs, fn_h5, param_chop):
# a function to create VAE dataset from a list of syllable spectrogram files
# first count the number of total windows so h5 can be initated
    def countChopWin_v2(fn):
        # distinguish two cases: 
        # 1. for syllables longer than chop size, count how many spectrogram chops
        # 2. for syllables shorter than chop size, will embed it in the middle, count 1 
        with h5py.File(fn, 'r') as f:
            spec = f['spec'][:]
        # how many spectrogram chops?
        if spec.shape[1]>=param_chop['delta']:
            chop_start = np.arange(0, spec.shape[1]-param_chop['delta']+1, param_chop['hop'])
            return(chop_start.shape[0])
        else:
            return(0)
    with Parallel(n_jobs=num_cpu, verbose=0) as parallel:
        results = parallel(delayed(countChopWin_v2)(fn) for fn in spec_dirs)
    img_cnt = np.sum(results)
    print(f'total number of spec windows: {img_cnt}')
    # read in one spectrogram to determine the dimension
    with h5py.File(spec_dirs[0], 'r') as f:
        spec_input = f['spec'][:]
    # # then create the h5 file
    with h5py.File(fn_h5, 'w') as hf:
        hf.create_dataset('spec_data', shape=(img_cnt, spec_input.shape[0], param_chop['delta']))
        count = 0
        id_pd = pd.DataFrame()  # save the id of each spectrogram chop in a dataframe
        for n in tqdm.tqdm(range(len(spec_dirs))):
            fn = spec_dirs[n]
            with h5py.File(fn, 'r') as f:
                spec  = f['spec'][:]
            if spec.shape[1]>=param_chop['delta']:  # syllable is equal or longer than chop size: chop
                chop_start = np.arange(0, spec.shape[1]-param_chop['delta']+1, param_chop['hop'])
                for i in range(len(chop_start)):
                    spec_this = spec[:, chop_start[i]:(chop_start[i]+param_chop['delta'])]
                    # add to the large h5 file
                    hf['spec_data'][count,:,:] = spec_this
                    id_this = pd.DataFrame([{'fn_h5':fn, 'idx_start':chop_start[i], 'idx_end': chop_start[i]+param_chop['delta']}])
                    id_pd = pd.concat([id_pd, id_this], ignore_index=True)
                    count += 1

    return(count, id_pd)


def SpecFocalSyllableSave_v5(fn, fd_save_this, param_chop, focal_syl):
    # given a wav file and its segment time/label files in the same wav folder, calculate the spectrograms of focal syllables
    # save syllable spectrograms into a separate folder
    # load the wav signal and associated segmentations
    data, fs = librosa.load(fn, sr=None)
    segments = np.loadtxt(fn.replace('.wav', '.time.txt'), delimiter=',', dtype='int')
    labels = np.loadtxt(fn.replace('.wav', '.label.txt'), delimiter=',', dtype='str', encoding='utf-8-sig')
    # only calculate for the focal syllables
    focal_idx = [ii for ii in range(labels.shape[0]) if labels[ii] in focal_syl]
    # ignore syllables that are shorter than one spectrogram frame
    thres_p = param_chop['nperseg']
    # dump the spectrograms of syllables into pickle files in a separate folder
    # define a within-scope function to save spectrograms of syllables
    def saveSegSpec(m):
         # calculate spectrograms of the segment only
        if (segments[m,1]-segments[m,0])>=thres_p:
            spec_tmp, dt, f, t = ZZ_get_spec(segments[m,0], segments[m,1], data, fs, nperseg=param_chop['nperseg'], noverlap=param_chop['noverlap'], 
                                                       target_freqs=param_chop['target_freqs'], target_time_dim=param_chop['target_time_dim'],
                                                       bool_clip=param_chop['bool_clip'], clim=param_chop['clim'])
            tmp_data = {}
            # save spectrogram windows as h5 file
            fn_h5 = os.path.join(fd_save_this, f'{os.path.basename(fn)}_{m:08d}_{labels[m]}.h5')
            # print(fn_h5)
            with h5py.File(fn_h5, 'w') as hf:
                hf.create_dataset('spec', data=spec_tmp)
                hf.create_dataset('fn_wav', data=fn)
                hf.create_dataset('start_idx', data=segments[m,0])
                hf.create_dataset('end_idx', data=segments[m,1])
                hf.create_dataset('fs', data=fs)
                hf.create_dataset('f', data=f)
                hf.create_dataset('t', data=t)
                hf.create_dataset('dt', data=dt)
                hf.create_dataset('filename', data=fn)
                hf.create_dataset('label', data=list(labels[m]))
                hf.create_dataset('syl_id', data=m)
    with Parallel(n_jobs=num_cpu, verbose=0) as parallel:
        parallel(delayed(saveSegSpec)(m) for m in focal_idx)
    return fd_save_this


def SpecFocalSyllableSave_v6(fn, fd_save_this, param_chop, focal_syl, suffix_time='.time.txt', suffix_label='.label.txt'):
    # given a wav file and its segment time/label files in the same wav folder, calculate the spectrograms of focal syllables
    # save syllable spectrograms into a separate folder
    # load the wav signal and associated segmentations
    data, fs = librosa.load(fn, sr=None)
    segments = np.loadtxt(fn.replace('.wav', suffix_time), delimiter=',', dtype='int')
    labels = np.loadtxt(fn.replace('.wav', suffix_label), delimiter=',', dtype='str', encoding='utf-8-sig')
    # only calculate for the focal syllables
    focal_idx = [ii for ii in range(labels.shape[0]) if labels[ii] in focal_syl]
    # ignore syllables that are shorter than one spectrogram frame
    thres_p = param_chop['nperseg']
    # dump the spectrograms of syllables into pickle files in a separate folder
    # define a within-scope function to save spectrograms of syllables
    def saveSegSpec(m):
         # calculate spectrograms of the segment only
        if (segments[m,1]-segments[m,0])>=thres_p:
            spec_tmp, dt, f, t = ZZ_get_spec(segments[m,0], segments[m,1], data, fs, nperseg=param_chop['nperseg'], noverlap=param_chop['noverlap'], 
                                                       target_freqs=param_chop['target_freqs'], target_time_dim=param_chop['target_time_dim'],
                                                       bool_clip=param_chop['bool_clip'], clim=param_chop['clim'])
            tmp_data = {}
            # save spectrogram windows as h5 file
            fn_h5 = os.path.join(fd_save_this, f'{os.path.basename(fn)}_{m:08d}_{labels[m]}.h5')
            # print(fn_h5)
            with h5py.File(fn_h5, 'w') as hf:
                hf.create_dataset('spec', data=spec_tmp)
                hf.create_dataset('fn_wav', data=fn)
                hf.create_dataset('start_idx', data=segments[m,0])
                hf.create_dataset('end_idx', data=segments[m,1])
                hf.create_dataset('fs', data=fs)
                hf.create_dataset('f', data=f)
                hf.create_dataset('t', data=t)
                hf.create_dataset('dt', data=dt)
                hf.create_dataset('filename', data=fn)
                hf.create_dataset('label', data=list(labels[m]))
                hf.create_dataset('syl_id', data=m)
    with Parallel(n_jobs=num_cpu, verbose=0) as parallel:
        parallel(delayed(saveSegSpec)(m) for m in focal_idx)
    return fd_save_this


def ZZ_createDatasetFocal_v1(spec_dirs, fn_h5, param_chop, win_width):
    # a function to create VAE dataset from a list of syllable spectrogram files
    # embed 
    img_cnt = len(spec_dirs)
    print(f'total number of spec windows: {img_cnt}')
    # read in one spectrogram to determine the dimension
    with h5py.File(spec_dirs[0], 'r') as f:
        spec_input = f['spec'][:]
    # then create the h5 file
    with h5py.File(fn_h5, 'w') as hf:
        hf.create_dataset('spec_data', shape=(img_cnt, spec_input.shape[0], win_width))
        id_pd = pd.DataFrame()  # save the id of each spectrogram chop in a dataframe
        for n in tqdm.tqdm(range(len(spec_dirs))):
            fn = spec_dirs[n]
            with h5py.File(fn, 'r') as f:
                spec  = f['spec'][:]
            # embed the syllables in the center of the fixed window
            win = np.zeros((spec_input.shape[0], win_width), dtype=spec.dtype)
            # calculate an offset on the start index
            off_idx = math.floor((win.shape[1]-spec.shape[1])/2)
            if off_idx>=0:  # if syllable is shorter than the fixed window, embed in the middle
                win[:, off_idx:(off_idx+spec.shape[1])] = spec
            else:  # if syllable is longer, take the central part of the syllable
                win = spec[:, (-off_idx):(-off_idx+win.shape[1])]
            # add to the large h5 file
            hf['spec_data'][n,:,:] = win
            id_this = pd.DataFrame([{'fn_h5':fn, 'idx_start':off_idx, 'idx_end': off_idx+spec.shape[1]}])
            id_pd = pd.concat([id_pd, id_this], ignore_index=True)
    count = len(spec_dirs)
    return(count, id_pd)


def SpecFocalSyllableSaveDbase_v1(fn_dbase, fd_save_this, param_chop, focal_syl):
    # load the dbase
    dbase = loadmat(fn_dbase, simplify_cells=True)
    dbase = dbase['dbase']
    # find all the entries that have focal_syl labels
    have_focal = [ii for ii in range(len(dbase['SoundFiles'])) if (focal_syl in dbase['SegmentTitles'][ii])]
    def saveSegSpec2(si):
        # for a given sound file index, load the audio file, calculate spectrogram
        fn = os.path.join(dbase['SoundFiles'][si]['folder'], dbase['SoundFiles'][si]['name'])
        # determine if it's production or auditory
        bAud_idx = np.where(dbase['PropertyNames']=='bAud')[0][0]
        if dbase['Properties'][si, bAud_idx]>0:
            aud_ch = 'chan17'
        else:
            aud_ch = 'chan0'
        fn_aud = fn.replace('chan0', aud_ch)
        # load the audio data
        dataset = netCDF4.Dataset(fn_aud, mode='r')
        data = dataset.variables['data'][:] / 10  # nc to wav divided by 10
        fs = dbase['Fs']
        segments = dbase['SegmentTimes'][si]
        labels = dbase['SegmentTitles'][si]
        if type(labels)==str:  # deal with special case where only 1 syllable in the file
            labels = np.array([labels])
            segments = np.expand_dims(segments, axis=0)
        # ignore syllables that are shorter than one spectrogram frame
        thres_p = param_chop['nperseg']
        syl_idx = np.where(labels==focal_syl)[0]
        for m in syl_idx:
            if (segments[m,1]-segments[m,0])>=thres_p:
                spec_tmp, dt, f, t = ZZ_get_spec(segments[m,0], segments[m,1], data, fs, nperseg=param_chop['nperseg'], noverlap=param_chop['noverlap'], 
                                                           target_freqs=param_chop['target_freqs'], target_time_dim=param_chop['target_time_dim'],
                                                           bool_clip=param_chop['bool_clip'], clim=param_chop['clim'])
                tmp_data = {}
                # save spectrogram windows as h5 file
                fn_h5 = os.path.join(fd_save_this, f'{aud_ch}_{os.path.basename(fn)}_{m:08d}_{labels[m]}.h5')
                # print(fn_h5)
                with h5py.File(fn_h5, 'w') as hf:
                    hf.create_dataset('spec', data=spec_tmp)
                    hf.create_dataset('fn_wav', data=fn)
                    hf.create_dataset('start_idx', data=segments[m,0])
                    hf.create_dataset('end_idx', data=segments[m,1])
                    hf.create_dataset('fs', data=fs)
                    hf.create_dataset('f', data=f)
                    hf.create_dataset('t', data=t)
                    hf.create_dataset('dt', data=dt)
                    hf.create_dataset('filename', data=fn)
                    hf.create_dataset('label', data=list(labels[m]))
                    hf.create_dataset('syl_id', data=m)
    with Parallel(n_jobs=num_cpu, verbose=0) as parallel:
        parallel(delayed(saveSegSpec2)(si) for si in have_focal)
    return fd_save_this


def SpecWinSaveDbase_v1(fn_dbase, fd_save_this, param_chop):
    # given a dbase that has WhisperSeg labels, calculate then save spectrograms
    # load the dbase
    dbase = loadmat(fn_dbase, simplify_cells=True)
    dbase = dbase['dbase']
    # find all file entries that have whisperSeg labels
    have_focal = [ii for ii in range(len(dbase['SoundFiles'])) if (len(dbase['SegmentTitles'][ii])>0)]
    # ignore syllables that are shorter than one spectrogram window
    thres_p = int(param_chop['nperseg']+(param_chop['delta']-1)*(param_chop['nperseg']-param_chop['noverlap']))
    def saveSegSpec3(si):
        # for a given sound file index, load the audio file, calculate spectrogram
        fn = os.path.join(dbase['SoundFiles'][si]['folder'], dbase['SoundFiles'][si]['name'])
        # determine if it's production or auditory
        bAud_idx = np.where(dbase['PropertyNames']=='bAud')[0][0]
        if dbase['Properties'][si, bAud_idx]>0:
            aud_ch = 'chan17'
        else:
            aud_ch = 'chan0'
        fn_aud = fn.replace('chan0', aud_ch)
        # load the audio data
        dataset = netCDF4.Dataset(fn_aud, mode='r')
        data = dataset.variables['data'][:] / 10  # nc to wav divided by 10
        fs = dbase['Fs']
        segments = dbase['SegmentTimes'][si]
        labels = dbase['SegmentTitles'][si]
        if type(labels)==str:  # deal with special case where only 1 syllable in the file
            labels = np.array([labels])
            segments = np.expand_dims(segments, axis=0)
        for m in range(segments.shape[0]):
            if (segments[m,1]-segments[m,0])>=thres_p:
                spec_tmp, dt, f, t = ZZ_get_spec(segments[m,0], segments[m,1], data, fs, nperseg=param_chop['nperseg'], noverlap=param_chop['noverlap'], 
                                                           target_freqs=param_chop['target_freqs'], target_time_dim=param_chop['target_time_dim'],
                                                           bool_clip=param_chop['bool_clip'], clim=param_chop['clim'])
                tmp_data = {}
                # save spectrogram windows as h5 file
                fn_h5 = os.path.join(fd_save_this, f'{os.path.basename(fn_aud)}_{m:08d}_{labels[m]}.h5')
                # print(fn_h5)
                with h5py.File(fn_h5, 'w') as hf:
                    hf.create_dataset('spec', data=spec_tmp)
                    hf.create_dataset('fn_wav', data=fn_aud)
                    hf.create_dataset('start_idx', data=segments[m,0])
                    hf.create_dataset('end_idx', data=segments[m,1])
                    hf.create_dataset('fs', data=fs)
                    hf.create_dataset('f', data=f)
                    hf.create_dataset('t', data=t)
                    hf.create_dataset('dt', data=dt)
                    hf.create_dataset('filename', data=fn_aud)
                    hf.create_dataset('label', data=list(labels[m]))
                    hf.create_dataset('syl_id', data=m)
                    
    with Parallel(n_jobs=num_cpu, verbose=0) as parallel:
        parallel(delayed(saveSegSpec3)(si) for si in have_focal)
        
        
def SpecWinSaveDbase_v2(fn_dbase, fd_save_this, param_chop, aud_select):
    # given a dbase that has WhisperSeg labels, calculate then save spectrograms
    # load the dbase
    dbase = loadmat(fn_dbase, simplify_cells=True)
    dbase = dbase['dbase']
    # find all file entries that have whisperSeg labels
    have_focal = [ii for ii in range(len(dbase['SoundFiles'])) if (len(dbase['SegmentTitles'][ii])>0)]
    # ignore syllables that are shorter than one spectrogram window
    thres_p = int(param_chop['nperseg']+(param_chop['delta']-1)*(param_chop['nperseg']-param_chop['noverlap']))
    def saveSegSpec3(si):
        # for a given sound file index, load the audio file, calculate spectrogram
        fn = os.path.join(dbase['SoundFiles'][si]['folder'], dbase['SoundFiles'][si]['name'])
        # determine if it's production or auditory
        bAud_idx = np.where(dbase['PropertyNames']=='bAud')[0][0]
        if dbase['Properties'][si, bAud_idx]>0:
            aud_ch = 'chan17'
        else:
            aud_ch = 'chan0'
        if aud_ch==aud_select:
            fn_aud = fn.replace('chan0', aud_ch)
            # load the audio data
            dataset = netCDF4.Dataset(fn_aud, mode='r')
            data = dataset.variables['data'][:] / 10  # nc to wav divided by 10
            fs = dbase['Fs']
            segments = dbase['SegmentTimes'][si]
            labels = dbase['SegmentTitles'][si]
            if type(labels)==str:  # deal with special case where only 1 syllable in the file
                labels = np.array([labels])
                segments = np.expand_dims(segments, axis=0)
            for m in range(segments.shape[0]):
                if (segments[m,1]-segments[m,0])>=thres_p:
                    spec_tmp, dt, f, t = ZZ_get_spec(segments[m,0], segments[m,1], data, fs, nperseg=param_chop['nperseg'], noverlap=param_chop['noverlap'], 
                                                               target_freqs=param_chop['target_freqs'], target_time_dim=param_chop['target_time_dim'],
                                                               bool_clip=param_chop['bool_clip'], clim=param_chop['clim'])
                    tmp_data = {}
                    # save spectrogram windows as h5 file
                    fn_h5 = os.path.join(fd_save_this, f'{os.path.basename(fn_aud)}_{m:08d}_{labels[m]}.h5')
                    # print(fn_h5)
                    with h5py.File(fn_h5, 'w') as hf:
                        hf.create_dataset('spec', data=spec_tmp)
                        hf.create_dataset('fn_wav', data=fn_aud)
                        hf.create_dataset('start_idx', data=segments[m,0])
                        hf.create_dataset('end_idx', data=segments[m,1])
                        hf.create_dataset('fs', data=fs)
                        hf.create_dataset('f', data=f)
                        hf.create_dataset('t', data=t)
                        hf.create_dataset('dt', data=dt)
                        hf.create_dataset('filename', data=fn_aud)
                        hf.create_dataset('label', data=list(labels[m]))
                        hf.create_dataset('syl_id', data=m)
                    
    with Parallel(n_jobs=num_cpu, verbose=0) as parallel:
        parallel(delayed(saveSegSpec3)(si) for si in have_focal)
        
        

def SpecWinSaveDbase_v3(fn_dbase, fd_save_this, param_chop, aud_select, focal_syl):
    # given a dbase that has WhisperSeg labels, calculate then save spectrograms of syllables according to parameters in param_chop 
    # differ from v2: can input a list of syllable labels focal_syl and corresponding list of what audio channels to use
    # e.g. focal_syl=['v', 'q', 'p'], aud_select=['bSorted', 'chan0', 'chan17'];  
    # aud_select='bSorted' means use the dbase Property field to determine whether to use 'chan0' or 'chan17'
    # load the dbase
    dbase = loadmat(fn_dbase, simplify_cells=True)
    dbase = dbase['dbase']
    fs = dbase['Fs']
    bAud_idx = np.where(dbase['PropertyNames']=='bAud')[0][0]
    bSorted_idx = np.where(dbase['PropertyNames']=='bSorted')[0][0]
    bProd_idx = np.where(dbase['PropertyNames']=='bProd')[0][0]
    # find all file entries that have whisperSeg labels
    have_focal = [ii for ii in range(len(dbase['SoundFiles'])) if (len(dbase['SegmentTitles'][ii])>0)]
    # ignore syllables that are shorter than one spectrogram frame
    thres_p = param_chop['nperseg']
    
    # given a sound file index, calculate spectrograms for focal syllables
    def saveSegSpec4(si):
        # si = 169
        # check if this sound file has the focal syllables
        segments = dbase['SegmentTimes'][si]
        labels = dbase['SegmentTitles'][si]
        # convert to 2d array if only element
        if len(labels)==1:
            segments = segments[np.newaxis,:]
        # find the index of the focal syllable
        syl_idx = [ii for ii in range(len(labels)) if labels[ii] in focal_syl]
        if syl_idx:
            # load the audio file from both focal and partner bird
            fn = os.path.join(dbase['SoundFiles'][si]['folder'], dbase['SoundFiles'][si]['name'])
            dataset = netCDF4.Dataset(fn, mode='r')
            if not (hasattr(dataset, 'variables') and 'data' in dataset.variables.keys()):
                return
            d_focal = dataset.variables['data'][:] / 10  # nc to wav divided by 10
            fn_p = fn.replace('chan0.nc', 'chan17.nc')
            dataset = netCDF4.Dataset(fn_p, mode='r')
            if not (hasattr(dataset, 'variables') and 'data' in dataset.variables.keys()):
                return
            d_partner = dataset.variables['data'][:] / 10  # nc to wav divided by 10
            # loop through each focal syllable, determine what audio file to use, then calculate spectrogram
            for m in syl_idx:
                if (segments[m,1]-segments[m,0])>=thres_p:            
                    syl_this = labels[m]
                    find_idx = focal_syl.index(syl_this)
                    aud_this = aud_select[find_idx]
                    aud_ch = ''
                    if aud_this=='bSorted':
                        # determine what audio file to use, only extract if it's sorted
                        if dbase['Properties'][si, bSorted_idx]>0 and dbase['Properties'][si, bAud_idx]>0:
                            aud_ch = 'chan17'
                        if dbase['Properties'][si, bSorted_idx]>0 and dbase['Properties'][si, bProd_idx]>0:
                            aud_ch = 'chan0'
                    else:
                        aud_ch = aud_this
                    if aud_ch:
                        if aud_ch=='chan17':
                            data = d_partner
                        else:
                            data = d_focal
                        spec_tmp, dt, f, t = ZZ_get_spec(segments[m,0], segments[m,1], data, fs, nperseg=param_chop['nperseg'], noverlap=param_chop['noverlap'], 
                                                                               target_freqs=param_chop['target_freqs'], target_time_dim=param_chop['target_time_dim'],
                                                                               bool_clip=param_chop['bool_clip'], clim=param_chop['clim'])
                        tmp_data = {}
                        # save spectrogram windows as h5 file
                        fd_save_now = os.path.join(fd_save_this, aud_ch)
                        fn_h5 = os.path.join(fd_save_now, f'{os.path.basename(fn)}_{m:08d}_{labels[m]}.h5')
                        # print(fn_h5)
                        with h5py.File(fn_h5, 'w') as hf:
                            hf.create_dataset('spec', data=spec_tmp)
                            hf.create_dataset('fn_wav', data=fn)
                            hf.create_dataset('start_idx', data=segments[m,0])
                            hf.create_dataset('end_idx', data=segments[m,1])
                            hf.create_dataset('fs', data=fs)
                            hf.create_dataset('f', data=f)
                            hf.create_dataset('t', data=t)
                            hf.create_dataset('dt', data=dt)
                            hf.create_dataset('filename', data=fn)
                            hf.create_dataset('label', data=list(labels[m]))
                            hf.create_dataset('syl_id', data=m)
    with Parallel(n_jobs=num_cpu, verbose=0) as parallel:
        parallel(delayed(saveSegSpec4)(si) for si in have_focal)