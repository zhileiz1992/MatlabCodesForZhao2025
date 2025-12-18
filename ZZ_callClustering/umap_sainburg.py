# this script contains functions used by Tim Sainburg to calculate spectrogram and run UMAP
# adapted from https://github.com/timsainb/avgn_paper
# Zhilei, 05/29/2025

import librosa, os, sys
import librosa.filters
import numpy as np
from scipy import signal
from scipy.signal import butter, lfilter
from librosa import mel_frequencies
from PIL import Image
import matplotlib.pyplot as plt
from tqdm.autonotebook import tqdm
import matplotlib.collections as mcoll
import matplotlib.path as mpath
from matplotlib import collections as mc
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib import gridspec
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import cKDTree
from matplotlib import lines
import matplotlib.colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd

# init hyperparameters
class HParams(object):
    """ Hparams was removed from tf 2.0alpha so this is a placeholder
    """
    def __init__(self, **kwargs):
        self.set_defaults()
        self.__dict__.update(kwargs)

    def set_defaults(self):
        self.win_length_ms = 5
        self.hop_length_ms = 1
        self.n_fft = 1024
        self.fs = 20000
        self.ref_level_db = 20
        self.min_level_db = -60
        self.preemphasis = 0.97
        self.mel = True
        self.num_mel_bins = 64
        self.mel_lower_edge_hertz = 200
        self.mel_upper_edge_hertz = 15000
        self.power = 1.5  # for spectral inversion
        self.griffin_lim_iters = 50
        self.butter_lowcut = 500
        self.butter_highcut = 15000
        self.reduce_noise = False
        self.noise_reduce_kwargs = {}
        self.mask_spec = False
        self.mask_spec_kwargs = {"spec_thresh": 0.9, "offset": 1e-10}
        self.nex = -1
        self.n_jobs = -1
        self.verbosity = 1


# make spectrograms
def _stft(y, fs, hparams):
    return librosa.stft(
        y=y,
        n_fft=hparams.n_fft,
        hop_length=int(hparams.hop_length_ms / 1000 * fs),
        win_length=int(hparams.win_length_ms / 1000 * fs),
    )

def _amp_to_db(x):
    return 20 * np.log10(np.maximum(1e-5, x))

def _normalize(S, hparams):
    return np.clip((S - hparams.min_level_db) / -hparams.min_level_db, 0, 1)

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

def spectrogram_nn(y, fs, hparams):
    D = _stft(preemphasis(y, hparams), fs, hparams)
    S = _amp_to_db(np.abs(D)) - hparams.ref_level_db
    return S

def preemphasis(x, hparams):
    return signal.lfilter([1, -hparams.preemphasis], [1], x)

def spectrogram(y, fs, hparams):
    return _normalize(spectrogram_nn(y, fs, hparams), hparams)

def norm(x):
    return (x - np.min(x)) / (np.max(x) - np.min(x))

def int16_to_float32(data):
    """ Converts from uint16 wav to float32 wav
    """
    if np.max(np.abs(data)) > 32768:
        raise ValueError("Data has values above 32768")
    return (data / 32768.0).astype("float32")


def make_spec(
    syll_wav,
    fs,
    hparams,
    mel_matrix=None,
    use_mel=True,
    norm_uint8=False,
):
    """
    """
    # convert to float
    if type(syll_wav[0]) == int:
        syll_wav = int16_to_float32(syll_wav)

    # create spec
    spec = spectrogram(syll_wav, fs, hparams)
    if use_mel:
        spec = np.dot(spec.T, mel_matrix).T
    if norm_uint8:
        spec = (norm(spec) * 255).astype("uint8")

    return spec


def prepare_mel_matrix(hparams, rate, return_numpy=True, GPU_backend=False):
    """ Create mel filter
    """
    # import tensorflow if needed
    if "tf" not in sys.modules:
        if not GPU_backend:
            os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"  # see issue #152
            os.environ["CUDA_VISIBLE_DEVICES"] = ""
        import tensorflow as tf

    # create a filter to convolve with the spectrogram
    mel_matrix = tf.signal.linear_to_mel_weight_matrix(
        num_mel_bins=hparams.num_mel_bins,
        num_spectrogram_bins=int(hparams.n_fft / 2) + 1,
        sample_rate=rate,
        lower_edge_hertz=hparams.mel_lower_edge_hertz,
        upper_edge_hertz=hparams.mel_upper_edge_hertz,
        dtype=tf.dtypes.float32,
        name=None,
    )

    # gets the center frequencies of mel bands
    mel_f = mel_frequencies(
        n_mels=hparams.num_mel_bins + 2,
        fmin=hparams.mel_lower_edge_hertz,
        fmax=hparams.mel_upper_edge_hertz,
    )

    # Slaney-style mel is scaled to be approx constant energy per channel (from librosa)
    enorm = tf.dtypes.cast(
        tf.expand_dims(
            tf.constant(
                2.0
                / (mel_f[2 : hparams.num_mel_bins + 2] - mel_f[: hparams.num_mel_bins])
            ),
            0,
        ),
        tf.float32,
    )

    mel_matrix = tf.multiply(mel_matrix, enorm)
    mel_matrix = tf.divide(mel_matrix, tf.reduce_sum(mel_matrix, axis=0))
    if return_numpy:
        return mel_matrix.numpy()
    else:
        return mel_matrix


def mask_spec(spec, spec_thresh=0.9, offset=1e-10):
    """ mask threshold a spectrogram to be above some % of the maximum power
    """
    mask = spec >= (spec.max(axis=0, keepdims=1) * spec_thresh + offset)
    return spec * mask


def log_resize_spec(spec, scaling_factor=10):
    resize_shape = [int(np.log(np.shape(spec)[1]) * scaling_factor), np.shape(spec)[0]]
    resize_spec = np.array(Image.fromarray(spec).resize(resize_shape, Image.Resampling.LANCZOS))
    return resize_spec


def pad_spectrogram(spectrogram, pad_length):
    """ Pads a spectrogram to being a certain length
    """
    excess_needed = pad_length - np.shape(spectrogram)[1]
    pad_left = np.floor(float(excess_needed) / 2).astype("int")
    pad_right = np.ceil(float(excess_needed) / 2).astype("int")
    return np.pad(
        spectrogram, [(0, 0), (pad_left, pad_right)], "constant", constant_values=0
    )

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


def flatten_spectrograms(specs):
    return np.reshape(specs, (np.shape(specs)[0], np.prod(np.shape(specs)[1:])))


## plotting functions
def scatter_spec(
    z,
    specs,
    column_size=10,
    pal_color="hls",
    matshow_kwargs={"cmap": plt.cm.Greys},
    scatter_kwargs={"alpha": 0.5, "s": 1},
    line_kwargs={"lw": 1, "ls": "dashed", "alpha": 1},
    color_points=False,
    figsize=(10, 10),
    range_pad=0.1,
    x_range=None,
    y_range=None,
    enlarge_points=0,
    draw_lines=True,
    n_subset=-1,
    ax=None,
    show_scatter=True,
    border_line_width=1,
    img_origin="lower",
):
    """
    """
    n_columns = column_size * 4 - 4
    pal = sns.color_palette(pal_color, n_colors=n_columns)

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(column_size, column_size)

    if x_range is None and y_range is None:
        xmin, xmax = np.sort(np.vstack(z)[:, 0])[
            np.array([int(len(z) * 0.01), int(len(z) * 0.99)])
        ]
        ymin, ymax = np.sort(np.vstack(z)[:, 1])[
            np.array([int(len(z) * 0.01), int(len(z) * 0.99)])
        ]
        # xmin, ymin = np.min(z, axis=0)
        # xmax, ymax = np.max(z, axis=0)
        xmin -= (xmax - xmin) * range_pad
        xmax += (xmax - xmin) * range_pad
        ymin -= (ymax - ymin) * range_pad
        ymax += (ymax - ymin) * range_pad
    else:
        xmin, xmax = x_range
        ymin, ymax = y_range

    x_block = (xmax - xmin) / column_size
    y_block = (ymax - ymin) / column_size

    # ignore segments outside of range
    z = np.array(z)
    mask = np.array(
        [(z[:, 0] > xmin) & (z[:, 1] > ymin) & (z[:, 0] < xmax) & (z[:, 1] < ymax)]
    )[0]

    if "labels" in scatter_kwargs:
        scatter_kwargs["labels"] = np.array(scatter_kwargs["labels"])[mask]
    specs = np.array(specs)[mask]
    z = z[mask]

    # prepare the main axis
    main_ax = fig.add_subplot(gs[1 : column_size - 1, 1 : column_size - 1])
    # main_ax.scatter(z[:, 0], z[:, 1], **scatter_kwargs)
    if show_scatter:
        scatter_projections(projection=z, ax=main_ax, fig=fig, **scatter_kwargs)

    # loop through example columns
    axs = {}
    for column in range(n_columns):
        # get example column location
        if column < column_size:
            row = 0
            col = column

        elif (column >= column_size) & (column < (column_size * 2) - 1):
            row = column - column_size + 1
            col = column_size - 1

        elif (column >= ((column_size * 2) - 1)) & (column < (column_size * 3 - 2)):
            row = column_size - 1
            col = column_size - 3 - (column - column_size * 2)
        elif column >= column_size * 3 - 3:
            row = n_columns - column
            col = 0

        axs[column] = {"ax": fig.add_subplot(gs[row, col]), "col": col, "row": row}
        # label subplot
        """axs[column]["ax"].text(
            x=0.5,
            y=0.5,
            s=column,
            horizontalalignment="center",
            verticalalignment="center",
            transform=axs[column]["ax"].transAxes,
        )"""

        # sample a point in z based upon the row and column
        xpos = xmin + x_block * col + x_block / 2
        ypos = ymax - y_block * row - y_block / 2
        # main_ax.text(x=xpos, y=ypos, s=column, color=pal[column])

        axs[column]["xpos"] = xpos
        axs[column]["ypos"] = ypos

    main_ax.set_xlim([xmin, xmax])
    main_ax.set_ylim([ymin, ymax])

    # create a voronoi diagram over the x and y pos points
    points = [[axs[i]["xpos"], axs[i]["ypos"]] for i in axs.keys()]

    voronoi_kdtree = cKDTree(points)
    vor = Voronoi(points)

    # plot voronoi
    # voronoi_plot_2d(vor, ax = main_ax);

    # find where each point lies in the voronoi diagram
    z = z[:n_subset]
    point_dist, point_regions = voronoi_kdtree.query(list(z))

    lines_list = []
    # loop through regions and select a point
    for key in axs.keys():
        # sample a point in (or near) voronoi region
        nearest_points = np.argsort(np.abs(point_regions - key))
        possible_points = np.where(point_regions == point_regions[nearest_points][0])[0]
        chosen_point = np.random.choice(a=possible_points, size=1)[0]
        point_regions[chosen_point] = 1e4
        # plot point
        if enlarge_points > 0:
            if color_points:
                color = pal[key]
            else:
                color = "k"
            main_ax.scatter(
                [z[chosen_point, 0]],
                [z[chosen_point, 1]],
                color=color,
                s=enlarge_points,
            )
        # draw spec
        axs[key]["ax"].matshow(
            specs[chosen_point],
            origin=img_origin,
            interpolation="none",
            aspect="auto",
            **matshow_kwargs,
        )

        axs[key]["ax"].set_xticks([])
        axs[key]["ax"].set_yticks([])
        if color_points:
            plt.setp(axs[key]["ax"].spines.values(), color=pal[key])

        for i in axs[key]["ax"].spines.values():
            i.set_linewidth(border_line_width)

        # draw a line between point and image
        if draw_lines:
            mytrans = (
                axs[key]["ax"].transAxes + axs[key]["ax"].figure.transFigure.inverted()
            )

            line_end_pos = [0.5, 0.5]

            if axs[key]["row"] == 0:
                line_end_pos[1] = 0
            if axs[key]["row"] == column_size - 1:
                line_end_pos[1] = 1

            if axs[key]["col"] == 0:
                line_end_pos[0] = 1
            if axs[key]["col"] == column_size - 1:
                line_end_pos[0] = 0

            infig_position = mytrans.transform(line_end_pos)

            xpos, ypos = main_ax.transLimits.transform(
                (z[chosen_point, 0], z[chosen_point, 1])
            )

            mytrans2 = main_ax.transAxes + main_ax.figure.transFigure.inverted()
            infig_position_start = mytrans2.transform([xpos, ypos])

            color = pal[key] if color_points else "k"
            lines_list.append(
                lines.Line2D(
                    [infig_position_start[0], infig_position[0]],
                    [infig_position_start[1], infig_position[1]],
                    color=color,
                    transform=fig.transFigure,
                    **line_kwargs,
                )
            )
    if draw_lines:
        for l in lines_list:
            fig.lines.append(l)

    gs.update(wspace=0, hspace=0)
    # gs.update(wspace=0.5, hspace=0.5)

    fig = plt.gcf()

    if ax is not None:
        buf = io.BytesIO()
        plt.savefig(buf, dpi=300, bbox_inches="tight", pad_inches=0)
        buf.seek(0)
        im = Image.open(buf)
        ax.imshow(im)
        plt.close(fig)

    return fig, axs, main_ax, [xmin, xmax, ymin, ymax]



def scatter_projections(
    syllables=None,
    projection=None,
    labels=None,
    ax=None,
    figsize=(10, 10),
    alpha=0.1,
    s=1,
    color="k",
    color_palette="tab20",
    categorical_labels=True,
    show_legend=True,
    tick_pos="bottom",
    tick_size=16,
    cbar_orientation="vertical",
    log_x=False,
    log_y=False,
    grey_unlabelled=True,
    fig=None,
    colornorm=False,
    rasterized=True,
    equalize_axes=True,
    print_lab_dict=False,  # prints color scheme
):
    """ creates a scatterplot of syllables using some projection
    """
    if projection is None:
        if syllables is None:
            raise ValueError("Either syllables or projections must by passed")

        syllables_flattened = np.reshape(
            syllables, (np.shape(syllables)[0], np.prod(np.shape(syllables)[1:]))
        )

        # if no projection is passed, assume umap
        fit = umap.UMAP(min_dist=0.25, verbose=True)
        u_all = fit.fit_transform(syllables_flattened)

    # color labels
    if labels is not None:
        if categorical_labels:
            if (color_palette == "tab20") & (len(np.unique(labels)) < 20):
                pal = sns.color_palette(color_palette, n_colors=20)
                pal = np.array(pal)[
                    np.linspace(0, 19, len(np.unique(labels))).astype("int")
                ]
                # print(pal)
            else:
                pal = sns.color_palette(color_palette, n_colors=len(np.unique(labels)))
            lab_dict = {lab: pal[i] for i, lab in enumerate(np.unique(labels))}
            if grey_unlabelled:
                if -1 in lab_dict.keys():
                    lab_dict[-1] = [0.95, 0.95, 0.95, 1.0]
                if print_lab_dict:
                    print(lab_dict)
            colors = np.array([lab_dict[i] for i in labels])
    else:
        colors = color

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

        # plot
    if colornorm:
        norm = norm = matplotlib.colors.LogNorm()
    else:
        norm = None
    if categorical_labels or labels is None:
        ax.scatter(
            projection[:, 0],
            projection[:, 1],
            rasterized=rasterized,
            alpha=alpha,
            s=s,
            color=colors,
            norm=norm,
        )

    else:
        cmin = np.quantile(labels, 0.01)
        cmax = np.quantile(labels, 0.99)
        sct = ax.scatter(
            projection[:, 0],
            projection[:, 1],
            vmin=cmin,
            vmax=cmax,
            cmap=plt.get_cmap(color_palette),
            rasterized=rasterized,
            alpha=alpha,
            s=s,
            c=labels,
        )

    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")

    if labels is not None:
        if categorical_labels == True:
            legend_elements = [
                Line2D([0], [0], marker="o", color=value, label=key)
                for key, value in lab_dict.items()
            ]
        if show_legend:
            if not categorical_labels:
                if cbar_orientation == "horizontal":
                    axins1 = inset_axes(
                        ax,
                        width="50%",  # width = 50% of parent_bbox width
                        height="5%",  # height : 5%
                        loc="upper left",
                    )
                    # cbar = fig.colorbar(sct, cax=axins1, orientation=cbar_orientation

                else:
                    axins1 = inset_axes(
                        ax,
                        width="5%",  # width = 50% of parent_bbox width
                        height="50%",  # height : 5%
                        loc="lower right",
                    )
                cbar = fig.colorbar(sct, cax=axins1, orientation=cbar_orientation)
                cbar.ax.tick_params(labelsize=tick_size)
                axins1.xaxis.set_ticks_position(tick_pos)
            else:
                ax.legend(handles=legend_elements)
    if equalize_axes:
        ax.axis("equal")
    return ax


## ZZ utility fuction
def ZZ_specFromWav_v1(fn, hparams, syl, mel_matrix): 
    # given wav file and selected syllables, output spectrograms and meta info
    specs = []
    info = pd.DataFrame()
    # check if annotation file exist
    fn_label = fn.replace('.wav','.label.txt')
    if os.path.exists(fn_label) and os.path.getsize(fn_label)>0:
        labels = np.genfromtxt(fn_label, dtype=str, delimiter=None, encoding='utf-8-sig')
        labels = np.atleast_1d(labels)
        # check if it has target syllable
        idx = [ii for ii in range(labels.shape[0]) if labels[ii] in syl]
        if len(idx)>0:
            fn_time = fn_label.replace('.label.txt', '.time.txt')
            seg = np.loadtxt(fn_time, delimiter=',', dtype=int)
            seg = np.atleast_2d(seg)
            # load wav file
            data, fs = librosa.load(fn, sr=None)
            # convert data to float32 if needed
            if np.issubdtype(type(data[0]), np.integer):
                data = int16_to_float32(data)
            # bandpass filter
            data = butter_bandpass_filter(data, hparams.butter_lowcut, hparams.butter_highcut, fs, order=5)
            # loop through target syllable
            for si in idx:
                x = data[seg[si,0]:seg[si,1]]
                # make spectrogram: with mel scaling
                spec = make_spec(x, fs, hparams, mel_matrix=mel_matrix, use_mel=hparams.mel, norm_uint8=False)
                # spec = umap_sainburg.make_spec(x, fs, hparams, mel_matrix=mel_matrix, use_mel=False, norm_uint8=False)
                specs.append(spec)
                # record the meta info of this syllable
                row = pd.DataFrame([{'fn_wav':fn, 's_idx':si, 'istart':seg[si,0], 'iend':seg[si,1], 'label':labels[si]}])
                info = pd.concat([info, row], ignore_index=True)
    return(specs, info)
