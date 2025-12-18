#!/usr/bin/env python
# coding: utf-8

# ## 0. Goal
# Given a dbase with WhisperSeg annotations, apply trained VAE/UMAP/HDBSCAN models to get labels for call subtypes

# In[1]:


import os, sys, importlib, librosa, glob, h5py, tqdm, pickle
from scipy.io import wavfile
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from joblib import Parallel, delayed
import random
import umap, hdbscan
from collections import Counter
import seaborn as sns
from matplotlib.colors import ListedColormap
from sklearn.metrics import silhouette_score
from torch.utils.data import Dataset, DataLoader
import torch
from sklearn.metrics import calinski_harabasz_score

plt.rcParams['pdf.fonttype'] = 42 


# In[2]:


# import my utility script
cluster_script_path = '/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_callClustering/'
sys.path.insert(1, cluster_script_path)
import hopkins2
import vae_goffinet, hopkins
importlib.reload(vae_goffinet)
importlib.reload(hopkins2)


# In[3]:


# create a custom colormap for spectrogram
jet = plt.get_cmap('jet', 255)
# Extract jet colors and prepend black at the beginning
jet_colors = jet(np.linspace(0, 1, 255))
custom_colors = np.vstack([[0, 0, 0, 1], jet_colors])  # Black for 0, then jet
custom_cmap = ListedColormap(custom_colors)


# In[ ]:





# ## 1. dbase -> spectrogram H5 datasets

# In[16]:


fn_dbase = FN_DBASE
fn_vae = FN_VAE
fn_umap = FN_UMAP
fn_hdbscan = FN_HDBSCAN
fd_temp = FD_TEMP
clim = CLIM
max_dur = MAX_DUR
syl = SYL

fd_save = fd_temp


# In[8]:


# define parameters for spectrograms
X_SHAPE = [128, 128]
p = {
    'get_spec': vae_goffinet.get_specZZ, # spectrogram maker
    'max_dur': 1e9, # maximum syllable duration
    'min_freq': 250, # minimum frequency
    'max_freq': 7500, # maximum frequency, default 7500
    'num_freq_bins': X_SHAPE[0], # hard-coded
    'num_time_bins': X_SHAPE[1], # hard-coded
    'nperseg': 256, # FFT
    'noverlap': 176, # FFT, determines window overlap when calculating spectrograms
    'spec_min_val': clim[0], # minimum log-spectrogram value
    'spec_max_val': clim[1], # maximum log-spectrogram value
    'fs': 20000, # audio samplerate
    'mel': False, # frequency spacing, mel or linear
    'time_stretch': False, # stretch short syllables?
    'within_syll_normalize': False, # normalize spectrogram values on a # spectrogram-by-spectrogram basis
    'max_num_syllables': None, # maximum number of syllables per directory
    'sylls_per_file': 20, # syllable per file, not used
    'real_preprocess_params': ('min_freq', 'max_freq', 'spec_min_val', 'spec_max_val'), # tunable parameters
    'int_preprocess_params': ('nperseg','noverlap'), # tunable parameters
    'binary_preprocess_params': ('time_stretch', 'mel', 'within_syll_normalize'), # tunable parameters
    'window_length': 0.032,  # for continuous sliding: size of the spectrogram window, unit is sec
    'hop_length': 0.004  # for continuous sliding: how much to hop for successive spectrogram windows, unit is sec
}


# In[15]:


# save spectrograms as .h5 files, along with meta info file
fn_h5 = os.path.join(fd_save, f'applyModel.{syl}.h5')
fn_info = os.path.join(fd_save, f'applyModel.{syl}.info.csv')
print(fn_h5)
print(fn_info)


# In[17]:


## grab all wav files
fns_wav = sorted(glob.glob(os.path.join(fd_temp, '*.wav')))
fns_label = sorted(glob.glob(os.path.join(fd_temp, '*.label.txt')))
print(f'Total number of wav files: {len(fns_wav)}')


# In[18]:


## calculate spectrogram
print('Calculating spectrograms...')
with Parallel(n_jobs=48, verbose=5) as parallel:
    res = parallel(delayed(vae_goffinet.ZZ_specFromWavGoffinet_v1)(fn, p, syl, max_dur) for fn in fns_wav)


# In[19]:


# flatten the result
temp = [aa[0] for aa in res]
specs = [arr for sublist in temp if sublist for arr in sublist]
df_list = [aa[1] for aa in res]
info = pd.concat([df for df in df_list if not df.empty], ignore_index=True)
print(len(specs), info.shape)


# In[20]:


## save results
# save padded spectrograms as h5 file 
spec_win_all = np.stack(specs, axis=0)
print(spec_win_all.shape)
with h5py.File(fn_h5, 'w') as f:
    f.create_dataset('spec_win_all', data=spec_win_all)

# save meta info as well
info.to_csv(fn_info)


# In[22]:


# plot some example spectrograms
# nrow = 3
# ncol = 8
# random.seed(1118)
# idx_rd = random.sample(range(len(specs)), nrow*ncol)
# fig, axes = plt.subplots(nrow, ncol, figsize=(14, nrow*2.5))
# for ii in range(len(idx_rd)):
#     plot_i = ii//ncol
#     plot_j = ii%ncol
#     ax = axes[plot_i][plot_j]
#     ax.imshow(specs[idx_rd[ii]], aspect='auto', cmap=custom_cmap, vmin=0, vmax=1, origin='lower')
#     ax.set_title(info.loc[idx_rd[ii], 'label'])
#     # set y tick labels
#     query_freqs = [1000, 3000, 5000, 7000]
#     target_freqs = info.loc[idx_rd[ii], 'spec_f']
#     indices = np.arange(len(target_freqs))
#     # Interpolate: given a value, find where it lies in the index space
#     interp_indices = np.interp(query_freqs, target_freqs, indices)
#     ax.set_yticks(interp_indices)
#     ax.set_yticklabels(query_freqs)
# plt.tight_layout()
# plt.show()
# fn_fig = os.path.join(fd_save, f'applyModel.exampleSpec.pdf')
# fig.savefig(fn_fig)


# In[ ]:





# ## 2. Apply trained VAE model

# In[23]:


model = vae_goffinet.VAE(save_dir=fd_save)
model.load_state(fn_vae)


# In[24]:


## Obtain latent representation
train_data = vae_goffinet.SpecDataset(fn_h5)
train_dataloader = DataLoader(train_data, batch_size=64, shuffle=False, num_workers=4)  # set shuffle to false to match the order in id_pd

# loop through dataloader, obtain model latent space
latent_m = np.zeros((info.shape[0], 32))
latent_d = np.zeros((info.shape[0], 32))
recon = np.zeros((info.shape[0], X_SHAPE[0], X_SHAPE[1]))
model.eval()
count = 0
for i, data in tqdm.tqdm(enumerate(train_dataloader)):
    data = data.to('cuda:0')
    with torch.no_grad():
        _, _, rec, mu, d = model.forwardZZ(data, return_latent_rec=True)
        a = rec.shape[0]
        latent_m[count:(count+a),:] = mu
        latent_d[count:(count+a),:] = d
        recon[count:(count+a),:,:] = rec
        count += a

# save the latent representations
fn_latentM = os.path.join(fd_save, 'latentM.csv')
np.savetxt(fn_latentM, latent_m, delimiter=',')
fn_latentD = os.path.join(fd_save, 'latentD.csv')
np.savetxt(fn_latentD, latent_d, delimiter=',')
print(latent_m.shape)


# In[26]:


## check recontruction accuracy (optional)
# plot some random samples
# fig, ax = plt.subplots(2, 10, figsize=[12,4])
# random.seed(1118)
# random_i = random.sample(list(range(recon.shape[0])), 10)
# with h5py.File(fn_h5, 'r') as file:
#     for ii in range(10):
#         spec = train_data[random_i[ii],:,:].numpy()
#         ax[0][ii].imshow(np.flip(spec, 0), aspect='auto', vmin=0, vmax=1, cmap=custom_cmap)
#         spec = recon[random_i[ii],:,:]
#         ax[1][ii].imshow(np.flip(spec, 0), aspect='auto', vmin=0, vmax=1, cmap=custom_cmap)
#         if ii>0:
#             ax[0][ii].axis('off')
#             ax[1][ii].axis('off')
# plt.tight_layout()
# # save fig
# fn_fig = os.path.join(fd_save, 'reconstructed_spectrogram.pdf')
# fig.savefig(fn_fig)


# In[ ]:





# ## 3. Apply trained UMAP model

# In[28]:


with open(fn_umap, 'rb') as f:
    UMAP_model = pickle.load(f)


# In[29]:


res= UMAP_model.transform(latent_m)


# In[30]:


# add latent_m to the embedding data frame
embed = info.copy()
for ii in range(latent_m.shape[1]):
    embed[f'vae{ii}'] = latent_m[:,ii]
# add UMAP embedding to the dataframe
for jj in range(res.shape[1]):
    embed[f'umap{jj+1}'] = res[:,jj]


# In[ ]:





# ## 4. Apply the HDBSCAN model

# In[38]:


with open(fn_hdbscan, 'rb') as f:
    cluster = pickle.load(f)


# In[39]:


test_labels, prob = hdbscan.approximate_predict(cluster, res)


# In[42]:


# add to the embed data frame
embed['hdbscan_cluster'] = test_labels + 1    # no clustering is label 0
embed['hdbscan_prob'] = prob

# save embedding and clustering 
fn_embed = os.path.join(fd_save, f'applyModel.{syl}.embedding.csv')
embed.to_csv(fn_embed)


# In[43]:


fn_embed


# In[ ]:




