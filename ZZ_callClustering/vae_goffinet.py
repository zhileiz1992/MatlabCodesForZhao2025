# Try the original VAE architecture in Goffinet...Pearson 2021
# Directly adapt from the codes in AVA package: https://github.com/pearsonlab/autoencoded-vocal-analysis

import numpy as np
import os
import torch
from torch.distributions import LowRankMultivariateNormal
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam
from torch.utils.data import Dataset, DataLoader
import warnings
from scipy.signal import stft
from scipy.interpolate import RegularGridInterpolator
import h5py
import pandas as pd
from scipy.io import wavfile
from scipy.signal import butter, lfilter
import random
from skimage import transform
import umap


X_SHAPE = (128,128)
"""Processed spectrogram shape: ``[freq_bins, time_bins]``"""
X_DIM = np.prod(X_SHAPE)
"""Processed spectrogram dimension: ``freq_bins * time_bins``"""


class VAE(nn.Module):
    """Variational Autoencoder class for single-channel images.

    Attributes
    ----------
    save_dir : str, optional
        Directory where the model is saved. Defaults to ``''``.
    lr : float, optional
        Model learning rate. Defaults to ``1e-3``.
    z_dim : int, optional
        Latent dimension. Defaults to ``32``.
    model_precision : float, optional
        Precision of the observation model. Defaults to ``10.0``.
    device_name : {'cpu', 'cuda', 'auto'}, optional
        Name of device to train the model on. When ``'auto'`` is passed,
        ``'cuda'`` is chosen if ``torch.cuda.is_available()``, otherwise
        ``'cpu'`` is chosen. Defaults to ``'auto'``.

    Notes
    -----
    The model is trained to maximize the standard ELBO objective:

    .. math:: \mathcal{L} = \mathbb{E}_{q(z|x)} log p(x,z) + \mathbb{H}[q(z|x)]

    where :math:`p(x,z) = p(z)p(x|z)` and :math:`\mathbb{H}` is differential
    entropy. The prior :math:`p(z)` is a unit spherical normal distribution. The
    conditional distribution :math:`p(x|z)` is set as a spherical normal
    distribution to prevent overfitting. The variational distribution,
    :math:`q(z|x)` is an approximately rank-1 multivariate normal distribution.
    Here, :math:`q(z|x)` and :math:`p(x|z)` are parameterized by neural
    networks. Gradients are passed through stochastic layers via the
    reparameterization trick, implemented by the PyTorch `rsample` method.

    The dimensions of the network are hard-coded for use with 128 x 128
    spectrograms. Although a desired latent dimension can be passed to
    `__init__`, the dimensions of the network limit the practical range of
    values roughly 8 to 64 dimensions. Fiddling with the image dimensions will
    require updating the parameters of the layers defined in `_build_network`.
    """

    def __init__(self, save_dir='', lr=1e-3, z_dim=32, model_precision=10.0,
        device_name="auto"):
        """Construct a VAE.

        Parameters
        ----------
        save_dir : str, optional
            Directory where the model is saved. Defaults to the current working
            directory.
        lr : float, optional
            Learning rate of the ADAM optimizer. Defaults to 1e-3.
        z_dim : int, optional
            Dimension of the latent space. Defaults to 32.
        model_precision : float, optional
            Precision of the noise model, p(x|z) = N(mu(z), \Lambda) where
            \Lambda = model_precision * I. Defaults to 10.0.
        device_name: str, optional
            Name of device to train the model on. Valid options are ["cpu",
            "cuda", "auto"]. "auto" will choose "cuda" if it is available.
            Defaults to "auto".

        Note
        ----
        - The model is built before it's parameters can be loaded from a file.
            This means `self.z_dim` must match `z_dim` of the model being
            loaded.
        """
        super(VAE, self).__init__()
        self.save_dir = save_dir
        self.lr = lr
        self.z_dim = z_dim
        self.model_precision = model_precision
        assert device_name != "cuda" or torch.cuda.is_available()
        if device_name == "auto":
            device_name = "cuda" if torch.cuda.is_available() else "cpu"
        self.device = torch.device(device_name)
        if self.save_dir != '' and not os.path.exists(self.save_dir):
            os.makedirs(self.save_dir)
        self._build_network()
        self.optimizer = Adam(self.parameters(), lr=self.lr)
        self.epoch = 0
        self.loss = {'train':{}, 'test':{}}
        self.to(self.device)


    def _build_network(self):
        """Define all the network layers."""
        # Encoder
        self.conv1 = nn.Conv2d(1, 8, 3,1,padding=1)
        self.conv2 = nn.Conv2d(8, 8, 3,2,padding=1)
        self.conv3 = nn.Conv2d(8, 16,3,1,padding=1)
        self.conv4 = nn.Conv2d(16,16,3,2,padding=1)
        self.conv5 = nn.Conv2d(16,24,3,1,padding=1)
        self.conv6 = nn.Conv2d(24,24,3,2,padding=1)
        self.conv7 = nn.Conv2d(24,32,3,1,padding=1)
        self.bn1 = nn.BatchNorm2d(1)
        self.bn2 = nn.BatchNorm2d(8)
        self.bn3 = nn.BatchNorm2d(8)
        self.bn4 = nn.BatchNorm2d(16)
        self.bn5 = nn.BatchNorm2d(16)
        self.bn6 = nn.BatchNorm2d(24)
        self.bn7 = nn.BatchNorm2d(24)
        self.fc1 = nn.Linear(8192,1024)
        self.fc2 = nn.Linear(1024,256)
        self.fc31 = nn.Linear(256,64)
        self.fc32 = nn.Linear(256,64)
        self.fc33 = nn.Linear(256,64)
        self.fc41 = nn.Linear(64,self.z_dim)
        self.fc42 = nn.Linear(64,self.z_dim)
        self.fc43 = nn.Linear(64,self.z_dim)
        # Decoder
        self.fc5 = nn.Linear(self.z_dim,64)
        self.fc6 = nn.Linear(64,256)
        self.fc7 = nn.Linear(256,1024)
        self.fc8 = nn.Linear(1024,8192)
        self.convt1 = nn.ConvTranspose2d(32,24,3,1,padding=1)
        self.convt2 = nn.ConvTranspose2d(24,24,3,2,padding=1,output_padding=1)
        self.convt3 = nn.ConvTranspose2d(24,16,3,1,padding=1)
        self.convt4 = nn.ConvTranspose2d(16,16,3,2,padding=1,output_padding=1)
        self.convt5 = nn.ConvTranspose2d(16,8,3,1,padding=1)
        self.convt6 = nn.ConvTranspose2d(8,8,3,2,padding=1,output_padding=1)
        self.convt7 = nn.ConvTranspose2d(8,1,3,1,padding=1)
        self.bn8 = nn.BatchNorm2d(32)
        self.bn9 = nn.BatchNorm2d(24)
        self.bn10 = nn.BatchNorm2d(24)
        self.bn11 = nn.BatchNorm2d(16)
        self.bn12 = nn.BatchNorm2d(16)
        self.bn13 = nn.BatchNorm2d(8)
        self.bn14 = nn.BatchNorm2d(8)


    def _get_layers(self):
        """Return a dictionary mapping names to network layers."""
        return {'fc1':self.fc1, 'fc2':self.fc2, 'fc31':self.fc31,
                'fc32':self.fc32, 'fc33':self.fc33, 'fc41':self.fc41,
                'fc42':self.fc42, 'fc43':self.fc43, 'fc5':self.fc5,
                'fc6':self.fc6, 'fc7':self.fc7, 'fc8':self.fc8, 'bn1':self.bn1,
                'bn2':self.bn2, 'bn3':self.bn3, 'bn4':self.bn4, 'bn5':self.bn5,
                'bn6':self.bn6, 'bn7':self.bn7, 'bn8':self.bn8, 'bn9':self.bn9,
                'bn10':self.bn10, 'bn11':self.bn11, 'bn12':self.bn12,
                'bn13':self.bn13, 'bn14':self.bn14, 'conv1':self.conv1,
                'conv2':self.conv2, 'conv3':self.conv3, 'conv4':self.conv4,
                'conv5':self.conv5, 'conv6':self.conv6, 'conv7':self.conv7,
                'convt1':self.convt1, 'convt2':self.convt2,
                'convt3':self.convt3, 'convt4':self.convt4,
                'convt5':self.convt5, 'convt6':self.convt6,
                'convt7':self.convt7}


    def encode(self, x):
        """
        Compute :math:`q(z|x)`.

        .. math:: q(z|x) = \mathcal{N}(\mu, \Sigma)
        .. math:: \Sigma = u u^{T} + \mathtt{diag}(d)

        where :math:`\mu`, :math:`u`, and :math:`d` are deterministic functions
        of `x` and :math:`\Sigma` denotes a covariance matrix.

        Parameters
        ----------
        x : torch.Tensor
            The input images, with shape: ``[batch_size, height=128,
            width=128]``

        Returns
        -------
        mu : torch.Tensor
            Posterior mean, with shape ``[batch_size, self.z_dim]``
        u : torch.Tensor
            Posterior covariance factor, as defined above. Shape:
            ``[batch_size, self.z_dim]``
        d : torch.Tensor
            Posterior diagonal factor, as defined above. Shape:
            ``[batch_size, self.z_dim]``
        """
        x = x.unsqueeze(1)
        x = F.relu(self.conv1(self.bn1(x)))
        x = F.relu(self.conv2(self.bn2(x)))
        x = F.relu(self.conv3(self.bn3(x)))
        x = F.relu(self.conv4(self.bn4(x)))
        x = F.relu(self.conv5(self.bn5(x)))
        x = F.relu(self.conv6(self.bn6(x)))
        x = F.relu(self.conv7(self.bn7(x)))
        x = x.view(-1, 8192)
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        mu = F.relu(self.fc31(x))
        mu = self.fc41(mu)
        u = F.relu(self.fc32(x))
        u = self.fc42(u).unsqueeze(-1) # Last dimension is rank \Sigma = 1.
        d = F.relu(self.fc33(x))
        d = torch.exp(self.fc43(d)) # d must be positive.
        return mu, u, d


    def decode(self, z):
        """
        Compute :math:`p(x|z)`.

        .. math:: p(x|z) = \mathcal{N}(\mu, \Lambda)

        .. math:: \Lambda = \mathtt{model\_precision} \cdot I

        where :math:`\mu` is a deterministic function of `z`, :math:`\Lambda` is
        a precision matrix, and :math:`I` is the identity matrix.

        Parameters
        ----------
        z : torch.Tensor
            Batch of latent samples with shape ``[batch_size, self.z_dim]``

        Returns
        -------
        x : torch.Tensor
            Batch of means mu, described above. Shape: ``[batch_size,
            X_DIM=128*128]``
        """
        z = F.relu(self.fc5(z))
        z = F.relu(self.fc6(z))
        z = F.relu(self.fc7(z))
        z = F.relu(self.fc8(z))
        z = z.view(-1,32,16,16)
        z = F.relu(self.convt1(self.bn8(z)))
        z = F.relu(self.convt2(self.bn9(z)))
        z = F.relu(self.convt3(self.bn10(z)))
        z = F.relu(self.convt4(self.bn11(z)))
        z = F.relu(self.convt5(self.bn12(z)))
        z = F.relu(self.convt6(self.bn13(z)))
        z = self.convt7(self.bn14(z))
        return z.view(-1, X_DIM)


    def forward(self, x, return_latent_rec=False):
        """
        Send `x` round trip and compute a loss.

        In more detail: Given `x`, compute :math:`q(z|x)` and sample:
        :math:`\hat{z} \sim q(z|x)` . Then compute :math:`\log p(x|\hat{z})`,
        the log-likelihood of `x`, the input, given :math:`\hat{z}`, the latent
        sample. We will also need the likelihood of :math:`\hat{z}` under the
        model's prior: :math:`p(\hat{z})`, and the entropy of the latent
        conditional distribution, :math:`\mathbb{H}[q(z|x)]` . ELBO can then be
        estimated as:

        .. math:: \\frac{1}{N} \sum_{i=1}^N \mathbb{E}_{\hat{z} \sim q(z|x_i)}
            \log p(x_i,\hat{z}) + \mathbb{H}[q(z|x_i)]

        where :math:`N` denotes the number of samples from the data distribution
        and the expectation is estimated using a single latent sample,
        :math:`\hat{z}`. In practice, the outer expectation is estimated using
        minibatches.

        Parameters
        ----------
        x : torch.Tensor
            A batch of samples from the data distribution (spectrograms).
            Shape: ``[batch_size, height=128, width=128]``
        return_latent_rec : bool, optional
            Whether to return latent means and reconstructions. Defaults to
            ``False``.

        Returns
        -------
        loss : torch.Tensor
            Negative ELBO times the batch size. Shape: ``[]``
        latent : numpy.ndarray, if `return_latent_rec`
            Latent means. Shape: ``[batch_size, self.z_dim]``
        reconstructions : numpy.ndarray, if `return_latent_rec`
            Reconstructed means. Shape: ``[batch_size, height=128, width=128]``
        """
        mu, u, d = self.encode(x)
        latent_dist = LowRankMultivariateNormal(mu, u, d)
        z = latent_dist.rsample()
        x_rec = self.decode(z)
        # E_{q(z|x)} p(z)
        elbo = -0.5 * (torch.sum(torch.pow(z,2)) + self.z_dim * np.log(2*np.pi))
        # E_{q(z|x)} p(x|z)
        pxz_term = -0.5 * X_DIM * (np.log(2*np.pi/self.model_precision))
        l2s = torch.sum(torch.pow(x.view(x.shape[0],-1) - x_rec, 2), dim=1)
        pxz_term = pxz_term - 0.5 * self.model_precision * torch.sum(l2s)
        elbo = elbo + pxz_term
        # H[q(z|x)]
        elbo = elbo + torch.sum(latent_dist.entropy())
        if return_latent_rec:
            return -elbo, z.detach().cpu().numpy(), \
                x_rec.view(-1, X_SHAPE[0], X_SHAPE[1]).detach().cpu().numpy()
        return -elbo
    

    def forwardZZ(self, x, return_latent_rec=False):
        mu, u, d = self.encode(x)
        latent_dist = LowRankMultivariateNormal(mu, u, d)
        z = latent_dist.rsample()
        x_rec = self.decode(z)
        # E_{q(z|x)} p(z)
        elbo = -0.5 * (torch.sum(torch.pow(z,2)) + self.z_dim * np.log(2*np.pi))
        # E_{q(z|x)} p(x|z)
        pxz_term = -0.5 * X_DIM * (np.log(2*np.pi/self.model_precision))
        l2s = torch.sum(torch.pow(x.view(x.shape[0],-1) - x_rec, 2), dim=1)
        pxz_term = pxz_term - 0.5 * self.model_precision * torch.sum(l2s)
        elbo = elbo + pxz_term
        # H[q(z|x)]
        elbo = elbo + torch.sum(latent_dist.entropy())
        if return_latent_rec:
            return -elbo, z.detach().cpu().numpy(), \
                x_rec.view(-1, X_SHAPE[0], X_SHAPE[1]).detach().cpu().numpy(), \
                mu.detach().cpu().numpy(), d.detach().cpu().numpy()
        return -elbo


    def train_epoch(self, train_loader):
        """
        Train the model for a single epoch.

        Parameters
        ----------
        train_loader : torch.utils.data.Dataloader
            ava.models.vae_dataset.SyllableDataset Dataloader for training set

        Returns
        -------
        elbo : float
            A biased estimate of the ELBO, estimated using samples from
            `train_loader`.
        """
        self.train()
        train_loss = 0.0
        for batch_idx, data in enumerate(train_loader):
            self.optimizer.zero_grad()
            data = data.to(self.device)
            loss = self.forward(data)
            train_loss += loss.item()
            loss.backward()
            self.optimizer.step()
        train_loss /= len(train_loader.dataset)
        print('Epoch: {} Average loss: {:.4f}'.format(self.epoch, \
                train_loss))
        self.epoch += 1
        return train_loss


    def test_epoch(self, test_loader):
        """
        Test the model on a held-out test set, return an ELBO estimate.

        Parameters
        ----------
        test_loader : torch.utils.data.Dataloader
            ava.models.vae_dataset.SyllableDataset Dataloader for test set

        Returns
        -------
        elbo : float
            An unbiased estimate of the ELBO, estimated using samples from
            `test_loader`.
        """
        self.eval()
        test_loss = 0.0
        with torch.no_grad():
            for i, data in enumerate(test_loader):
                data = data.to(self.device)
                loss = self.forward(data)
                test_loss += loss.item()
        test_loss /= len(test_loader.dataset)
        print('Test loss: {:.4f}'.format(test_loss))
        return test_loss


    def train_loop(self, loaders, epochs=100, test_freq=2, save_freq=10,
        vis_freq=1):
        """
        Train the model for multiple epochs, testing and saving along the way.

        Parameters
        ----------
        loaders : dictionary
            Dictionary mapping the keys ``'test'`` and ``'train'`` to respective
            torch.utils.data.Dataloader objects.
        epochs : int, optional
            Number of (possibly additional) epochs to train the model for.
            Defaults to ``100``.
        test_freq : int, optional
            Testing is performed every `test_freq` epochs. Defaults to ``2``.
        save_freq : int, optional
            The model is saved every `save_freq` epochs. Defaults to ``10``.
        vis_freq : int, optional
            Syllable reconstructions are plotted every `vis_freq` epochs.
            Defaults to ``1``.
        """
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
            # Save the model.
            if (save_freq is not None) and (epoch % save_freq == 0) and \
                    (epoch > 0):
                filename = "checkpoint_"+str(epoch).zfill(3)+'.tar'
                self.save_state(filename)
            # Plot reconstructions.
            if (vis_freq is not None) and (epoch % vis_freq == 0):
                self.visualize(loaders['test'])


    def save_state(self, filename):
        """Save all the model parameters to the given file."""
        layers = self._get_layers()
        state = {}
        for layer_name in layers:
            state[layer_name] = layers[layer_name].state_dict()
        state['optimizer_state'] = self.optimizer.state_dict()
        state['loss'] = self.loss
        state['z_dim'] = self.z_dim
        state['epoch'] = self.epoch
        state['lr'] = self.lr
        state['save_dir'] = self.save_dir
        filename = os.path.join(self.save_dir, filename)
        torch.save(state, filename)


    def load_state(self, filename):
        """
        Load all the model parameters from the given ``.tar`` file.

        The ``.tar`` file should be written by `self.save_state`.

        Parameters
        ----------
        filename : str
            File containing a model state.

        Note
        ----
        - `self.lr`, `self.save_dir`, and `self.z_dim` are not loaded.
        """
        checkpoint = torch.load(filename, map_location=self.device)
        assert checkpoint['z_dim'] == self.z_dim
        layers = self._get_layers()
        for layer_name in layers:
            layer = layers[layer_name]
            layer.load_state_dict(checkpoint[layer_name])
        self.optimizer.load_state_dict(checkpoint['optimizer_state'])
        self.loss = checkpoint['loss']
        self.epoch = checkpoint['epoch']


    def get_latent(self, loader):
        """
        Get latent means for all syllable in the given loader.

        Parameters
        ----------
        loader : torch.utils.data.Dataloader
            ava.models.vae_dataset.SyllableDataset Dataloader.

        Returns
        -------
        latent : numpy.ndarray
            Latent means. Shape: ``[len(loader.dataset), self.z_dim]``

        Note
        ----
        - Make sure your loader is not set to shuffle if you're going to match
          these with labels or other fields later.
        """
        latent = np.zeros((len(loader.dataset), self.z_dim))
        i = 0
        for data in loader:
            data = data.to(self.device)
            with torch.no_grad():
                mu, _, _ = self.encode(data)
            mu = mu.detach().cpu().numpy()
            latent[i:i+len(mu)] = mu
            i += len(mu)
        return latent


## Han's data loader
class SpecDataset(Dataset):
    '''
    Grab spectrogram slices from a h5 file
    
    filepath = r'\train_data.h5'
    tmp_h5 = h5py.File(filepath, 'r')
    group_key = list(tmp_h5.keys())[0]
    ds_all =  tmp_h5[group_key]
    train_data= SpecDataset(ds_all )    
    
    '''
    def __init__(self, path):
        self.file_path = path
        self.dataset = None
        
        with h5py.File(self.file_path, 'r') as file:
            self.groupkey = list(file.keys())[0]
            self.dataset_len = len(file[self.groupkey])
        
    def __getitem__(self, index):
        if self.dataset is None:
            self.dataset = h5py.File(self.file_path, 'r')[self.groupkey]
        return torch.Tensor(self.dataset[index])
    
    def __len__(self):
        return self.dataset_len


# Pearson's script on making spectrograms
EPSILON = 1e-12
def get_specZZ(t1, t2, audio, p, fs=32000, target_freqs=None, target_times=None, \
    fill_value=-1/EPSILON, max_dur=None, remove_dc_offset=True):
    """
    Norm, scale, threshold, stretch, and resize a Short Time Fourier Transform.

    Notes
    -----
    * ``fill_value`` necessary?
    * Look at all references and see what can be simplified.
    * Why is a flag returned?

    Parameters
    ----------
    t1 : float
        Onset time.
    t2 : float
        Offset time.
    audio : numpy.ndarray
        Raw audio.
    p : dict
        Parameters. Must include keys: ...
    fs : float
        Samplerate.
    target_freqs : numpy.ndarray or ``None``, optional
        Interpolated frequencies.
    target_times : numpy.ndarray or ``None``, optional
        Intepolated times.
    fill_value : float, optional
        Defaults to ``-1/EPSILON``.
    max_dur : float, optional
        Maximum duration. Defaults to ``None``.
    remove_dc_offset : bool, optional
        Whether to remove any DC offset from the audio. Defaults to ``True``.

    Returns
    -------
    spec : numpy.ndarray
        Spectrogram.
    flag : bool
        ``True``
    """
    if max_dur is None:
        max_dur = p['max_dur']
    if t2 - t1 > max_dur + 1e-4:
        message = "Found segment longer than max_dur: " + str(t2-t1) + \
                "s, max_dur = " + str(max_dur) + "s"
        # warnings.warn(message)
    s1, s2 = int(round(t1*fs)), int(round(t2*fs))
    assert s1 < s2, "s1: " + str(s1) + " s2: " + str(s2) + " t1: " + str(t1) + \
            " t2: " + str(t2)
    # Get a spectrogram and define the interpolation object.
    temp = min(len(audio),s2) - max(0,s1)
    if temp < p['nperseg'] or s2 <= 0 or s1 >= len(audio):
        return np.zeros((p['num_freq_bins'], p['num_time_bins'])), True
    else:
        temp_audio = audio[max(0,s1):min(len(audio),s2)]
        if remove_dc_offset:
            temp_audio = temp_audio - np.mean(temp_audio)
        f, t, spec = stft(temp_audio, fs=fs, nperseg=p['nperseg'], \
                noverlap=p['noverlap'])
    t += max(0,t1)
    spec = np.log(np.abs(spec) + EPSILON)
    # interp = interp2d(t, f, spec, copy=False, bounds_error=False, \
    # 	fill_value=fill_value)
    interp = RegularGridInterpolator((f, t), spec, bounds_error=False, fill_value=fill_value)
    # Define target frequencies.
    if target_freqs is None:
        if p['mel']:
            target_freqs = np.linspace(_mel(p['min_freq']), \
                    _mel(p['max_freq']), p['num_freq_bins'])
            target_freqs = _inv_mel(target_freqs)
        else:
            target_freqs = np.linspace(p['min_freq'], p['max_freq'], \
                    p['num_freq_bins'])
    # Define target times.
    if target_times is None:
        duration = t2 - t1
        if p['time_stretch']:
            duration = np.sqrt(duration * max_dur) # stretched duration
        shoulder = 0.5 * (max_dur - duration)
        target_times = np.linspace(t1-shoulder, t2+shoulder, p['num_time_bins'])
    # Then interpolate.
    # interp_spec = interp(target_times, target_freqs, assume_sorted=True)
    F_mesh, T_mesh = np.meshgrid(target_freqs, target_times, indexing='ij')  # F_mesh.shape = (m, n)
    # Stack the meshgrid into points of shape (m*n, 2)
    points = np.stack([F_mesh.ravel(), T_mesh.ravel()], axis=-1)
    # Query the interpolator
    interpolated_values = interp(points)  # shape (m*n,)
    # Reshape the result back to 2D: shape (m, n)
    interp_spec = interpolated_values.reshape(F_mesh.shape)
    spec = interp_spec
    # Normalize.
    spec -= p['spec_min_val']
    spec /= (p['spec_max_val'] - p['spec_min_val'])
    spec = np.clip(spec, 0.0, 1.0)
    # Within-syllable normalize.
    if p['within_syll_normalize']:
        spec -= np.quantile(spec, p['normalize_quantile'])
        spec[spec<0.0] = 0.0
        spec /= np.max(spec) + EPSILON
    return spec, True, target_freqs, target_times





EPSILON = 1e-12
def get_specZZnoTstretch(t1, t2, audio, p, fs=32000, target_freqs=None, target_times=None, \
    fill_value=-1/EPSILON, max_dur=None, remove_dc_offset=True):
    """
    Norm, scale, threshold, stretch, and resize a Short Time Fourier Transform.
    ZZ: no interpolation on the time axis 
    Notes
    -----
    * ``fill_value`` necessary?
    * Look at all references and see what can be simplified.
    * Why is a flag returned?

    Parameters
    ----------
    t1 : float
        Onset time.
    t2 : float
        Offset time.
    audio : numpy.ndarray
        Raw audio.
    p : dict
        Parameters. Must include keys: ...
    fs : float
        Samplerate.
    target_freqs : numpy.ndarray or ``None``, optional
        Interpolated frequencies.
    target_times : numpy.ndarray or ``None``, optional
        Intepolated times.
    fill_value : float, optional
        Defaults to ``-1/EPSILON``.
    max_dur : float, optional
        Maximum duration. Defaults to ``None``.
    remove_dc_offset : bool, optional
        Whether to remove any DC offset from the audio. Defaults to ``True``.

    Returns
    -------
    spec : numpy.ndarray
        Spectrogram.
    flag : bool
        ``True``
    """
    if max_dur is None:
        max_dur = p['max_dur']
    if t2 - t1 > max_dur + 1e-4:
        message = "Found segment longer than max_dur: " + str(t2-t1) + \
                "s, max_dur = " + str(max_dur) + "s"
        warnings.warn(message)
    s1, s2 = int(round(t1*fs)), int(round(t2*fs))
    assert s1 < s2, "s1: " + str(s1) + " s2: " + str(s2) + " t1: " + str(t1) + \
            " t2: " + str(t2)
    # Get a spectrogram and define the interpolation object.
    temp = min(len(audio),s2) - max(0,s1)
    if temp < p['nperseg'] or s2 <= 0 or s1 >= len(audio):
        return np.zeros((p['num_freq_bins'], p['num_time_bins'])), True
    else:
        temp_audio = audio[max(0,s1):min(len(audio),s2)]
        if remove_dc_offset:
            temp_audio = temp_audio - np.mean(temp_audio)
        f, t, spec = stft(temp_audio, fs=fs, nperseg=p['nperseg'], \
                noverlap=p['noverlap'])
    t += max(0,t1)
    spec = np.log(np.abs(spec) + EPSILON)
    # interp = interp2d(t, f, spec, copy=False, bounds_error=False, \
    # 	fill_value=fill_value)
    interp = RegularGridInterpolator((f, t), spec, bounds_error=False, fill_value=fill_value)
    # Define target frequencies.
    if target_freqs is None:
        if p['mel']:
            target_freqs = np.linspace(_mel(p['min_freq']), \
                    _mel(p['max_freq']), p['num_freq_bins'])
            target_freqs = _inv_mel(target_freqs)
        else:
            target_freqs = np.linspace(p['min_freq'], p['max_freq'], \
                    p['num_freq_bins'])
    # Define target times.
    if target_times is None:
        target_times = t
    # Then interpolate.
    # interp_spec = interp(target_times, target_freqs, assume_sorted=True)
    F_mesh, T_mesh = np.meshgrid(target_freqs, target_times, indexing='ij')  # F_mesh.shape = (m, n)
    # Stack the meshgrid into points of shape (m*n, 2)
    points = np.stack([F_mesh.ravel(), T_mesh.ravel()], axis=-1)
    # Query the interpolator
    interpolated_values = interp(points)  # shape (m*n,)
    # Reshape the result back to 2D: shape (m, n)
    interp_spec = interpolated_values.reshape(F_mesh.shape)
    spec = interp_spec
    # Normalize.
    spec -= p['spec_min_val']
    spec /= (p['spec_max_val'] - p['spec_min_val'])
    spec = np.clip(spec, 0.0, 1.0)
    # Within-syllable normalize.
    if p['within_syll_normalize']:
        spec -= np.quantile(spec, p['normalize_quantile'])
        spec[spec<0.0] = 0.0
        spec /= np.max(spec) + EPSILON
    return spec, True, target_freqs, target_times


def _mel(a):
    """https://en.wikipedia.org/wiki/Mel-frequency_cepstrum"""
    return 1127 * np.log(1 + a / 700)


def _inv_mel(a):
    """https://en.wikipedia.org/wiki/Mel-frequency_cepstrum"""
    return 700 * (np.exp(a / 1127) - 1)


def get_syll_specsZZ(onsets, offsets, audio_filename, p):
    """
    Return the spectrograms corresponding to `onsets` and `offsets`.

    Parameters
    ----------
    onsets : list of floats
        Syllable onsets.
    offsets : list of floats
        Syllable offsets.
    audio_filename : str
        Audio filename.
    p : dict
        A dictionary mapping preprocessing parameters to their values. NOTE: ADD
        REFERENCE HERE!

    Returns
    -------
    specs : list of {numpy.ndarray, None}
        Spectrograms.
    valid_syllables : list of int
        Indices of `specs` containing valid syllables.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=WavFileWarning)
        fs, audio = wavfile.read(audio_filename)
    assert p['nperseg'] % 2 == 0 and p['nperseg'] > 2
    if p['mel']:
        target_freqs = np.linspace( \
                _mel(p['min_freq']), _mel(p['max_freq']), p['num_freq_bins'])
        target_freqs = _inv_mel(target_freqs)
    else:
        target_freqs = np.linspace( \
                p['min_freq'], p['max_freq'], p['num_freq_bins'])
    specs, valid_syllables = [], []
    # For each syllable...
    for i, t1, t2 in zip(range(len(onsets)), onsets, offsets):
        spec, valid = get_specZZ(t1, t2, audio, p, fs, \
                target_freqs=target_freqs)
        if valid:
            valid_syllables.append(i)
            specs.append(spec)
    return specs, valid_syllables

def flatten_spectrograms(specs):
    return np.reshape(specs, (np.shape(specs)[0], np.prod(np.shape(specs)[1:])))


def ZZ_getDur_v1(fn_label, fs, syl):
    dur = []
    if os.path.getsize(fn_label)>0:
        labels = np.genfromtxt(fn_label, dtype=str, delimiter=None, encoding='utf-8-sig')
        labels = np.atleast_1d(labels)  # deal with only one segment in a file
        # check if it has target syllable
        idx = [ii for ii in range(labels.shape[0]) if labels[ii] in syl]
        if len(idx)>0:
            fn_time = fn_label.replace('.label.txt', '.time.txt')
            seg = np.loadtxt(fn_time, delimiter=',', dtype=int)
            seg = np.atleast_2d(seg)
            dur = list((seg[idx,1] - seg[idx,0]) / fs)
    return(dur)


def ZZ_specFromWavGoffinet_v1(fn, p, syl, max_dur): 
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
            fs, audio = wavfile.read(fn)
            # butterworth filter
            audio = butter_bandpass_filter(audio, p['min_freq'], p['max_freq'], fs, order=5)
            duration = len(audio)/fs
            for m in idx:
                # determine onset and offset of syllables
                onset = seg[m,0] / fs
                offset = seg[m,1] / fs
                if (seg[m,1]-seg[m,0])<=p['nperseg']:
                    continue
                spec, flag, target_freqs, target_times = get_specZZ(onset, offset, audio, p, fs=fs, target_times=None, remove_dc_offset=True, max_dur=max_dur)
                specs.append(spec)
                # record the meta info of this syllable
                row = pd.DataFrame([{'fn_wav':fn, 's_idx':m, 'istart':seg[m,0], 'iend':seg[m,1], 'label':labels[m], 'spec_f':target_freqs, 'spec_t':target_times}])
                info = pd.concat([info, row], ignore_index=True)
    return(specs, info)


def ZZ_specFromWavGoffinet_v2(fn, p, syl, max_dur, dur_cutoff): 
    # calculate spectrograms from wav file
    # differ from v1: remove syllables with duraction shorter than dur_cutoff
    # check if annotation file exist
    fn_label = fn.replace('.wav','.label.txt')
    pt_cutoff = p['fs'] * dur_cutoff  # convert duration cutoff to point cutoff
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
            # check if duration pass cutoff
            idx2 = [ii for ii in idx if (seg[ii,1]-seg[ii,0])>=pt_cutoff]
            if len(idx2)>0: 
                # load wav file
                fs, audio = wavfile.read(fn)
                # butterworth filter
                audio = butter_bandpass_filter(audio, p['min_freq'], p['max_freq'], fs, order=5)
                for m in idx2:
                    # determine onset and offset of syllables
                    onset = seg[m,0] / fs
                    offset = seg[m,1] / fs
                    if (seg[m,1]-seg[m,0])<=p['nperseg']:
                        continue
                    spec, flag, target_freqs, target_times = get_specZZ(onset, offset, audio, p, fs=fs, target_times=None, remove_dc_offset=True, max_dur=max_dur)
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


EPSILON = 1e-12
def shuffle_specZZ(t1, t2, audio, p, fs=32000, target_freqs=None, target_times=None, fill_value=-1/EPSILON, max_dur=None, remove_dc_offset=True, block_size=3):
    """
    Norm, scale, threshold, stretch, and resize a Short Time Fourier Transform.
    Then shuffle the spectrogram

    Notes
    -----
    * ``fill_value`` necessary?
    * Look at all references and see what can be simplified.
    * Why is a flag returned?

    Parameters
    ----------
    t1 : float
        Onset time.
    t2 : float
        Offset time.
    audio : numpy.ndarray
        Raw audio.
    p : dict
        Parameters. Must include keys: ...
    fs : float
        Samplerate.
    target_freqs : numpy.ndarray or ``None``, optional
        Interpolated frequencies.
    target_times : numpy.ndarray or ``None``, optional
        Intepolated times.
    fill_value : float, optional
        Defaults to ``-1/EPSILON``.
    max_dur : float, optional
        Maximum duration. Defaults to ``None``.
    remove_dc_offset : bool, optional
        Whether to remove any DC offset from the audio. Defaults to ``True``.

    Returns
    -------
    spec : numpy.ndarray
        Spectrogram.
    flag : bool
        ``True``
    """
    if max_dur is None:
        max_dur = p['max_dur']
    if t2 - t1 > max_dur + 1e-4:
        message = "Found segment longer than max_dur: " + str(t2-t1) + \
                "s, max_dur = " + str(max_dur) + "s"
        # warnings.warn(message)
    s1, s2 = int(round(t1*fs)), int(round(t2*fs))
    assert s1 < s2, "s1: " + str(s1) + " s2: " + str(s2) + " t1: " + str(t1) + \
            " t2: " + str(t2)
    # Get a spectrogram and define the interpolation object.
    temp = min(len(audio),s2) - max(0,s1)
    # if temp < p['nperseg'] or s2 <= 0 or s1 >= len(audio):
    #     return np.zeros((p['num_freq_bins'], p['num_time_bins'])), True
    # else:
    temp_audio = audio[max(0,s1):min(len(audio),s2)]
    if remove_dc_offset:
        temp_audio = temp_audio - np.mean(temp_audio)
    f, t, spec = stft(temp_audio, fs=fs, nperseg=p['nperseg'], \
            noverlap=p['noverlap'])
    t += max(0,t1)
    spec = np.log(np.abs(spec) + EPSILON)
    # shuffle spec 
    spec = shuffle_blocks(spec, block_size)
    # interp = interp2d(t, f, spec, copy=False, bounds_error=False, \
    # 	fill_value=fill_value)
    interp = RegularGridInterpolator((f, t), spec, bounds_error=False, fill_value=fill_value)
    # Define target frequencies.
    if target_freqs is None:
        if p['mel']:
            target_freqs = np.linspace(_mel(p['min_freq']), \
                    _mel(p['max_freq']), p['num_freq_bins'])
            target_freqs = _inv_mel(target_freqs)
        else:
            target_freqs = np.linspace(p['min_freq'], p['max_freq'], \
                    p['num_freq_bins'])
    # Define target times.
    if target_times is None:
        duration = t2 - t1
        if p['time_stretch']:
            duration = np.sqrt(duration * max_dur) # stretched duration
        shoulder = 0.5 * (max_dur - duration)
        target_times = np.linspace(t1-shoulder, t2+shoulder, p['num_time_bins'])
    # Then interpolate.
    # interp_spec = interp(target_times, target_freqs, assume_sorted=True)
    F_mesh, T_mesh = np.meshgrid(target_freqs, target_times, indexing='ij')  # F_mesh.shape = (m, n)
    # Stack the meshgrid into points of shape (m*n, 2)
    points = np.stack([F_mesh.ravel(), T_mesh.ravel()], axis=-1)
    # Query the interpolator
    interpolated_values = interp(points)  # shape (m*n,)
    # Reshape the result back to 2D: shape (m, n)
    interp_spec = interpolated_values.reshape(F_mesh.shape)
    spec = interp_spec
    # Normalize.
    spec -= p['spec_min_val']
    spec /= (p['spec_max_val'] - p['spec_min_val'])
    spec = np.clip(spec, 0.0, 1.0)
    # Within-syllable normalize.
    if p['within_syll_normalize']:
        spec -= np.quantile(spec, p['normalize_quantile'])
        spec[spec<0.0] = 0.0
        spec /= np.max(spec) + EPSILON
    return spec, True, target_freqs, target_times


def shuffle_blocks(spec, block_size):
    n_rows, n_cols = spec.shape
    block_starts = list(range(0, n_cols, block_size))
    
    # Split into blocks
    blocks = [spec[:, start:start + block_size] for start in block_starts]
    
    # Shuffle the blocks in-place
    np.random.shuffle(blocks)
    
    # Concatenate the shuffled blocks
    shuffled_spec = np.concatenate(blocks, axis=1)
    return shuffled_spec


def ZZ_ShuffleSpecFromWavGoffinet_v1(fn, p, syl, max_dur, block_size): 
    # calculate spectrograms from wav file
    # then perform random shuffle on the time axis
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
            fs, audio = wavfile.read(fn)
            # butterworth filter
            audio = butter_bandpass_filter(audio, p['min_freq'], p['max_freq'], fs, order=5)
            duration = len(audio)/fs
            for m in idx:
                # determine onset and offset of syllables
                onset = seg[m,0] / fs
                offset = seg[m,1] / fs
                if (seg[m,1]-seg[m,0])<=p['nperseg']:
                    continue
                spec, flag, target_freqs, target_times = shuffle_specZZ(onset, offset, audio, p, fs=fs, target_times=None, remove_dc_offset=True, max_dur=max_dur, block_size=block_size)
                specs.append(spec)
                # record the meta info of this syllable
                row = pd.DataFrame([{'fn_wav':fn, 's_idx':m, 'istart':seg[m,0], 'iend':seg[m,1], 'label':labels[m], 'spec_f':target_freqs, 'spec_t':target_times}])
                info = pd.concat([info, row], ignore_index=True)
    return(specs, info)


def ZZ_ShuffleSpecFromWavGoffinet_v2(fn, p, syl, max_dur, block_size): 
    # calculate spectrograms from wav file
    # then perform one roll shift
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
            fs, audio = wavfile.read(fn)
            # butterworth filter
            audio = butter_bandpass_filter(audio, p['min_freq'], p['max_freq'], fs, order=5)
            duration = len(audio)/fs
            for m in idx:
                # determine onset and offset of syllables
                onset = seg[m,0] / fs
                offset = seg[m,1] / fs
                if (seg[m,1]-seg[m,0])<=p['nperseg']:
                    continue
                spec, flag, target_freqs, target_times = roll_specZZ(onset, offset, audio, p, fs=fs, target_times=None, remove_dc_offset=True, max_dur=max_dur)
                specs.append(spec)
                # record the meta info of this syllable
                row = pd.DataFrame([{'fn_wav':fn, 's_idx':m, 'istart':seg[m,0], 'iend':seg[m,1], 'label':labels[m], 'spec_f':target_freqs, 'spec_t':target_times}])
                info = pd.concat([info, row], ignore_index=True)
    return(specs, info)


EPSILON = 1e-12
def roll_specZZ(t1, t2, audio, p, fs=32000, target_freqs=None, target_times=None, fill_value=-1/EPSILON, max_dur=None, remove_dc_offset=True):
    """
    Norm, scale, threshold, stretch, and resize a Short Time Fourier Transform.
    Then shuffle the spectrogram

    Notes
    -----
    * ``fill_value`` necessary?
    * Look at all references and see what can be simplified.
    * Why is a flag returned?

    Parameters
    ----------
    t1 : float
        Onset time.
    t2 : float
        Offset time.
    audio : numpy.ndarray
        Raw audio.
    p : dict
        Parameters. Must include keys: ...
    fs : float
        Samplerate.
    target_freqs : numpy.ndarray or ``None``, optional
        Interpolated frequencies.
    target_times : numpy.ndarray or ``None``, optional
        Intepolated times.
    fill_value : float, optional
        Defaults to ``-1/EPSILON``.
    max_dur : float, optional
        Maximum duration. Defaults to ``None``.
    remove_dc_offset : bool, optional
        Whether to remove any DC offset from the audio. Defaults to ``True``.

    Returns
    -------
    spec : numpy.ndarray
        Spectrogram.
    flag : bool
        ``True``
    """
    if max_dur is None:
        max_dur = p['max_dur']
    if t2 - t1 > max_dur + 1e-4:
        message = "Found segment longer than max_dur: " + str(t2-t1) + \
                "s, max_dur = " + str(max_dur) + "s"
        # warnings.warn(message)
    s1, s2 = int(round(t1*fs)), int(round(t2*fs))
    assert s1 < s2, "s1: " + str(s1) + " s2: " + str(s2) + " t1: " + str(t1) + \
            " t2: " + str(t2)
    # Get a spectrogram and define the interpolation object.
    temp = min(len(audio),s2) - max(0,s1)
    # if temp < p['nperseg'] or s2 <= 0 or s1 >= len(audio):
    #     return np.zeros((p['num_freq_bins'], p['num_time_bins'])), True
    # else:
    temp_audio = audio[max(0,s1):min(len(audio),s2)]
    if remove_dc_offset:
        temp_audio = temp_audio - np.mean(temp_audio)
    f, t, spec = stft(temp_audio, fs=fs, nperseg=p['nperseg'], \
            noverlap=p['noverlap'])
    t += max(0,t1)
    spec = np.log(np.abs(spec) + EPSILON)
    # roll shift spec by a random shift value
    st = random.sample(range(spec.shape[1]),1)[0]
    spec = np.roll(spec, st, axis=1)
    # interp = interp2d(t, f, spec, copy=False, bounds_error=False, \
    # 	fill_value=fill_value)
    interp = RegularGridInterpolator((f, t), spec, bounds_error=False, fill_value=fill_value)
    # Define target frequencies.
    if target_freqs is None:
        if p['mel']:
            target_freqs = np.linspace(_mel(p['min_freq']), \
                    _mel(p['max_freq']), p['num_freq_bins'])
            target_freqs = _inv_mel(target_freqs)
        else:
            target_freqs = np.linspace(p['min_freq'], p['max_freq'], \
                    p['num_freq_bins'])
    # Define target times.
    if target_times is None:
        duration = t2 - t1
        if p['time_stretch']:
            duration = np.sqrt(duration * max_dur) # stretched duration
        shoulder = 0.5 * (max_dur - duration)
        target_times = np.linspace(t1-shoulder, t2+shoulder, p['num_time_bins'])
    # Then interpolate.
    # interp_spec = interp(target_times, target_freqs, assume_sorted=True)
    F_mesh, T_mesh = np.meshgrid(target_freqs, target_times, indexing='ij')  # F_mesh.shape = (m, n)
    # Stack the meshgrid into points of shape (m*n, 2)
    points = np.stack([F_mesh.ravel(), T_mesh.ravel()], axis=-1)
    # Query the interpolator
    interpolated_values = interp(points)  # shape (m*n,)
    # Reshape the result back to 2D: shape (m, n)
    interp_spec = interpolated_values.reshape(F_mesh.shape)
    spec = interp_spec
    # Normalize.
    spec -= p['spec_min_val']
    spec /= (p['spec_max_val'] - p['spec_min_val'])
    spec = np.clip(spec, 0.0, 1.0)
    # Within-syllable normalize.
    if p['within_syll_normalize']:
        spec -= np.quantile(spec, p['normalize_quantile'])
        spec[spec<0.0] = 0.0
        spec /= np.max(spec) + EPSILON
    return spec, True, target_freqs, target_times


EPSILON = 1e-12
def get_specZZrowT(temp_audio, p, fs=32000, target_freqs=None, fill_value=-1/EPSILON, remove_dc_offset=True):
    """
    Norm, scale, threshold, a Short Time Fourier Transform. No operation on time axis, return with original range and resolution
    Make spectrograms for shotgun VAE analysis

    Notes
    -----
    * ``fill_value`` necessary?
    * Look at all references and see what can be simplified.
    * Why is a flag returned?

    Parameters
    ----------
    temp_audio : numpy.ndarray
        Raw audio.
    p : dict
        Parameters. Must include keys: ...
    fs : float
        Samplerate.
    target_freqs : numpy.ndarray or ``None``, optional
        Interpolated frequencies.
    fill_value : float, optional
        Defaults to ``-1/EPSILON``.
    remove_dc_offset : bool, optional
        Whether to remove any DC offset from the audio. Defaults to ``True``.

    Returns
    -------
    spec : numpy.ndarray
        Spectrogram.
    flag : bool
        ``True``
    """

    # calculate power spectrogram 
    if remove_dc_offset:
        temp_audio = temp_audio - np.mean(temp_audio)
    f, t, spec = stft(temp_audio, fs=fs, nperseg=p['nperseg'], \
            noverlap=p['noverlap'])

    spec = np.log(np.abs(spec) + EPSILON)

    # define interpolator for frequency axis
    interp = RegularGridInterpolator((f, t), spec, bounds_error=False, fill_value=fill_value)
    # Define target frequencies.
    if target_freqs is None:
        if p['mel']:
            target_freqs = np.linspace(_mel(p['min_freq']), \
                    _mel(p['max_freq']), p['num_freq_bins'])
            target_freqs = _inv_mel(target_freqs)
        else:
            target_freqs = np.linspace(p['min_freq'], p['max_freq'], \
                    p['num_freq_bins'])
    # Define target times.
    target_times = t
    # Then interpolate.
    # interp_spec = interp(target_times, target_freqs, assume_sorted=True)
    F_mesh, T_mesh = np.meshgrid(target_freqs, target_times, indexing='ij')  # F_mesh.shape = (m, n)
    # Stack the meshgrid into points of shape (m*n, 2)
    points = np.stack([F_mesh.ravel(), T_mesh.ravel()], axis=-1)
    # Query the interpolator
    interpolated_values = interp(points)  # shape (m*n,)
    # Reshape the result back to 2D: shape (m, n)
    interp_spec = interpolated_values.reshape(F_mesh.shape)
    spec = interp_spec
    # Normalize.
    spec -= p['spec_min_val']
    spec /= (p['spec_max_val'] - p['spec_min_val'])
    spec = np.clip(spec, 0.0, 1.0)
    # Within-syllable normalize.
    if p['within_syll_normalize']:
        spec -= np.quantile(spec, p['normalize_quantile'])
        spec[spec<0.0] = 0.0
        spec /= np.max(spec) + EPSILON
    return spec, True, target_freqs, target_times


def ZZ_specFromWavTraj_v1(fn, p, syl): 
    # calculate spectrograms from wav file for shot-gun VAE
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
        # find all the syllables that belong to the specified list 
        if syl:
            idx = [ii for ii in range(labels.shape[0]) if labels[ii] in syl]
        else:
            idx = [ii for ii in range(labels.shape[0])]
        if len(idx)>0:
            fn_time = fn_label.replace('.label.txt', '.time.txt')
            seg = np.loadtxt(fn_time, delimiter=',', dtype=int)
            seg = np.atleast_2d(seg)
            # load wav file
            fs, audio = wavfile.read(fn)
            # how many data points to pad
            pad_pt = int(fs*p['pad'])
            # butterworth filter
            audio = butter_bandpass_filter(audio, p['min_freq'], p['max_freq'], fs, order=5)
            duration = len(audio)/fs
            for m in idx:
                onset = seg[m,0]
                offset = seg[m,1]
                # skip if too short for FFT
                if (seg[m,1]-seg[m,0])<=p['nperseg']:
                    continue
                # pad the onset/offset, if too close to boundary or another, pad additional zeros, such that syllable starts at pad_pt
                # also calculate how many additional zeros to pad at the start and end
                zero_start = 0
                zero_end = 0
                i_start = onset - pad_pt
                if (m==0) and (i_start<0):
                    zero_start = -i_start
                    i_start = 0
                # if bleed into previous syllable
                if (m>0) and (i_start<=seg[m-1,1]):
                    zero_start = seg[m-1,1] - i_start
                    i_start = seg[m-1, 1]
                i_end = offset + pad_pt
                if (m==seg.shape[0]-1) and i_end>=audio.shape[0]:
                    zero_end = i_end - audio.shape[0]
                    i_end = audio.shape[0]
                # if bleed into next syllable
                if (m<seg.shape[0]-1) and (i_end>=seg[m+1,0]):
                    zero_end = i_end - seg[m+1,0]
                    i_end = seg[m+1,0]
                # print(m, onset, offset, i_start, i_end, zero_start, zero_end)
                temp_audio = audio[i_start:i_end]
                # pad additional zeros if needed
                temp_audio = np.concatenate([np.zeros(zero_start, temp_audio.dtype), temp_audio, np.zeros(zero_end, temp_audio.dtype)])
                # plt.plot(temp_audio)
                # plt.xlim(0, 2000)
                spec, flag, target_freqs, target_times = get_specZZrowT(temp_audio, p, fs=fs, target_freqs=None, remove_dc_offset=True)
    
                specs.append(spec)
                # record the meta info of this syllable
                row = pd.DataFrame([{'fn_wav':fn, 's_idx':m, 'istart':seg[m,0], 'iend':seg[m,1], 'label':labels[m], 'spec_f':target_freqs, 'spec_t':target_times, 
                                    'i_start':i_start, 'i_end':i_end, 'zero_start':zero_start, 'zero_end':zero_end, 'rel_ori': pad_pt}])
                info = pd.concat([info, row], ignore_index=True)
    return(specs, info)


def ZZ_sampleSpecWin_v1(fn_spec, ri, num_win_per, p, resize=True):
    # a function to sample spectrogram windows from a syllable spectrogram dataset fn_spec and syllable index ri
    # use for making training dataset for VAE
    # get the spectrogram dataset
    with h5py.File(fn_spec, 'r') as f:
        # Access the dataset
        spec = f[f'spec_{ri}'][:] 

    # sample starting index for the sliding window
    win_frame = p['win_frame']
    hop_frame = p['hop_frame']
    win_pad = p['win_pad']
    # where does the syllable actually start
    hop_ms = (p['nperseg']-p['noverlap'])/p['fs']*1000
    onset_frame = int(p['pad']*1000 / hop_ms)
    # calculate the range where starting index can be sampled
    r_start = max([0, onset_frame-win_pad])
    # the end needs to consider the pad difference and the width of the sliding window
    r_end = spec.shape[1] - r_start - win_frame
    # determine how many sliding windows can be sampled
    act_num = min([num_win_per, r_end-r_start])

    # get the sampled spec window data
    # default to return an empty list
    spec_wins = []
    info_wins = pd.DataFrame()
    if act_num>0:
        i_start = random.sample(range(r_start,r_end), act_num)
        for i_i, i_s in enumerate(i_start):
            win_this = spec[:, i_s:(i_s+win_frame)] 
            # resize if requested
            if resize:
                win_this = transform.resize(win_this, (p['num_freq_bins'], p['num_time_bins']), order=1, mode='edge', anti_aliasing=True)
            spec_wins.append(win_this)
            row = pd.DataFrame([{'ri':ri, 'i_i':i_i, 'i_s':i_s, 'i_e':i_s+win_frame}])
            info_wins = pd.concat([info_wins, row], ignore_index=True)
    return(spec_wins, info_wins)


def ZZ_sampleSpecWin_v2(fn_spec, ri, p, win_prop=1, resize=True):
    # a function to sample spectrogram windows from a syllable spectrogram dataset fn_spec and syllable index ri
    # the number of windows to sample is proportional to the syllable duration: win_prop number of windows per sliding window duration of the syllable
    # use for making training dataset for VAE
    # get the spectrogram dataset
    with h5py.File(fn_spec, 'r') as f:
        # Access the dataset
        spec = f[f'spec_{ri}'][:] 

    # sample starting index for the sliding window
    win_frame = p['win_frame']
    hop_frame = p['hop_frame']
    win_pad = p['win_pad']
    # how many frames were padded for the syllable spectrogram 
    spec_pad = int(p['pad']*p['fs'] / (p['nperseg']-p['noverlap']))
    # where does the syllable actually start
    hop_ms = (p['nperseg']-p['noverlap'])/p['fs']*1000
    onset_frame = int(p['pad']*1000 / hop_ms)
    # calculate the range where starting index can be sampled
    r_start = max([0, onset_frame-win_pad])
    # the end needs to consider the pad difference and the width of the sliding window
    r_end = spec.shape[1] - r_start - win_frame
    # determine how many sliding windows to be sampled: proportional to the syllable duration
    num_win_per = int((spec.shape[1] - 2*spec_pad) // win_frame * win_prop)
    # sample one window if syllable duration is shorter than sliding window duration
    num_win_per = max([1, num_win_per])
    act_num = min([num_win_per, r_end-r_start])

    # get the sampled spec window data
    # default to return an empty list
    spec_wins = []
    info_wins = pd.DataFrame()
    if act_num>0:
        i_start = random.sample(range(r_start,r_end), act_num)
        for i_i, i_s in enumerate(i_start):
            win_this = spec[:, i_s:(i_s+win_frame)] 
            # resize if requested
            if resize:
                win_this = transform.resize(win_this, (p['num_freq_bins'], p['num_time_bins']), order=1, mode='edge', anti_aliasing=True)
            spec_wins.append(win_this)
            row = pd.DataFrame([{'ri':ri, 'i_i':i_i, 'i_s':i_s, 'i_e':i_s+win_frame}])
            info_wins = pd.concat([info_wins, row], ignore_index=True)
    return(spec_wins, info_wins)


def ZZ_slideSylWin_v1(fn_spec, ri, p, resize=True):
    # given a row index, load the spectrogram data, then chop 
    # use for making dataset for getting trajectory of an entire syllable
    with h5py.File(fn_spec, 'r') as f:
        # Access the dataset
        spec = f[f'spec_{ri}'][:] 

    # calculate the start/end index for sliding, consider the different pads
    win_frame = p['win_frame']
    hop_frame = p['hop_frame']
    win_pad = p['win_pad']
    # where does the syllable actually start
    hop_ms = (p['nperseg']-p['noverlap'])/p['fs']*1000
    onset_frame = int(p['pad']*1000 / hop_ms)
    # calculate the range where starting index can be sampled
    r_start = max([0, onset_frame-win_pad])
    # the end needs to consider the pad difference and the width of the sliding window
    r_end = spec.shape[1] - r_start - win_frame
    
    # perform slide
    # default to return an empty list
    spec_wins = []
    info_wins = pd.DataFrame()
    start_all = list(range(r_start, r_end, hop_frame))
    if len(start_all)>0:
        for i_i, i_s in enumerate(start_all):
            win_this = spec[:, i_s:(i_s+win_frame)] 
            # resize if requested
            if resize:
                win_this = transform.resize(win_this, (p['num_freq_bins'], p['num_time_bins']), order=1, mode='edge', anti_aliasing=True)
            spec_wins.append(win_this)
            row = pd.DataFrame([{'ri':ri, 'i_i':i_i, 'i_s':i_s, 'i_e':i_s+win_frame}])
            info_wins = pd.concat([info_wins, row], ignore_index=True)
            
    return(spec_wins, info_wins)


def ZZ_runUMAP_v1(d, param_umap, random_state=1118, meta_info=None):
    # a script to run UMAP with defined parameters, return the model and result
    umap_model = umap.UMAP(n_neighbors=param_umap['n_neighbors'], n_components=param_umap['n_components'], min_dist=param_umap['min_dist'], 
                                                  metric=param_umap['metric'], random_state=random_state, verbose=True)
    res = umap_model.fit_transform(d)

    # combine result with meta_info dataframe if specified
    if meta_info is not None:
        embed = meta_info.copy()
        for ii in range(d.shape[1]):
            embed[f'vae{ii}'] = d[:,ii]
        for jj in range(res.shape[1]):
            embed[f'umap{jj+1}'] = res[:,jj]
        return(umap_model, embed)
    else:
        return(umap_model, res)