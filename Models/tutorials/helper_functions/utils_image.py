#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import torch
import cv2
from plenoptic.simulate import SteerablePyramidFreq
import plenoptic as po
import einops
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from itertools import product
from scipy.ndimage import gaussian_filter

def get_stats_within_eccentricity(img, pixpdeg, min_ecc_deg=4, max_ecc_deg=8):
    """
    Compute min, max, and mean of pixel values within a given eccentricity range.

    Parameters
    ----------
    img : 2D numpy array
        Image (e.g., one frame of stim_energy[i,j][ch,h])
    pixpdeg : float
        Pixels per degree (used to convert eccentricity in degrees to pixels)
    min_ecc_deg : float
        Minimum eccentricity (deg) to include
    max_ecc_deg : float
        Maximum eccentricity (deg) to include

    Returns
    -------
    stats : dict
        {'min': ..., 'max': ..., 'mean': ...}
    """

    # Convert eccentricity limits to pixels
    min_ecc_px = min_ecc_deg * pixpdeg
    max_ecc_px = max_ecc_deg * pixpdeg

    # Image dimensions and center
    ny, nx = img.shape
    cy, cx = ny / 2, nx / 2

    # Distance of each pixel from the center (in pixels)
    y, x = np.ogrid[:ny, :nx]
    r = np.sqrt((x - cx)**2 + (y - cy)**2)

    # Mask for pixels within the eccentricity annulus
    mask = (r >= min_ecc_px) & (r <= max_ecc_px)

    # Extract values and compute stats
    vals = img[mask]
    stats = {
        'min': float(vals.min()),
        'max': float(vals.max()),
        'mean': float(vals.mean())
    }

    return stats

def retrieve_dim(stim_energy):
    # Get all keys and shapes
    keys = list(stim_energy.keys())
    sample_array = next(iter(stim_energy.values()))
    n_set = len(set(k[0] for k in keys))
    n_stimuli = len(set(k[1] for k in keys))
    n_ori = sample_array.shape[1]  # second dimension of your arrays
    n_sf = sample_array.shape[0]
    
    return n_set, n_stimuli, n_ori, n_sf

def visualize_filter_magnitudes(stim_energy, representative_chIdx, scale=None, share_color_range=False, title=None, cbar=False, scale_factor=3, save_path=None):

    [n_set, n_stimuli, n_ori, _] = retrieve_dim(stim_energy)
        
    #fig, axes = plt.subplots(n_set, n_stimuli, figsize=(n_stimuli*scale_factor, n_set*scale_factor))

    
    if scale is not None:
        vmin = scale[0]
        vmax = scale[1]
        if share_color_range is False:
            print("scale given, so assuming shared scale requested: changing share_color_range to True")
        share_color_range = True
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        if share_color_range:
            print("no scale given, so computing shared scale: global min/max")
            vmin = min(image.min() for image in stim_energy.values())  # Smallest value across all images
            vmax = max(image.max() for image in stim_energy.values())
            norm = colors.Normalize(vmin=vmin, vmax=vmax)
        else:
            print("no scale given, no shared scale requested")
            vmin = None
            vmax = None
    
    for i in range(n_set):  # e.g., Cartesian, Polar
    
        fig, axes = plt.subplots(
            n_ori, n_stimuli,
            figsize=(n_stimuli * 2.5, n_ori * 2.5),
            squeeze=False
        )
        fig.suptitle(f"Set {i}", fontsize=14)
    
        for h in range(n_ori):  # orientation
            for j in range(n_stimuli):  # stimulus index
                ax = axes[h, j]
    
                img = stim_energy[i, j][representative_chIdx, h]  # select correct image
                ax.axis('off')
    
                # shared vs independent color scale
                if share_color_range:
                    im = ax.imshow(img, cmap='gray', norm=norm)
                else:
                    im = ax.imshow(img, cmap='gray')
    
                # optional colorbar
                if cbar:
                    cb = plt.colorbar(im, ax=ax, format="%.2e", fraction=0.046, pad=0.04)
                    cb.ax.tick_params(labelsize=8)
    
        fig.supxlabel('Input Stimuli', fontsize=14)
        fig.supylabel('Orientation Channel', fontsize=14)
        
        #plt.tight_layout()
    
        plt.tight_layout(rect=[0.01, 0.01, 0.98, 0.90])  # leaves bottom & left space
    
        plt.subplots_adjust(top=0.8, wspace=0.03, hspace=0.03)
        fig.suptitle(title, fontsize=14, y=0.97)
        plt.show()
    
    if save_path is not None:
        fig.savefig(save_path, bbox_inches='tight', pad_inches=0)

    return vmin, vmax
        

def plot_energy_per_channel(stim_energy, representative_chIdx, pixpdeg, min_ecc_deg=4, max_ecc_deg=8, yrange=None):

    ori_labels = ["90°", "45°", "0°", "135°"]
    set_labels = ["Cartesian Set", "Polar Set"]
    
    metric="mean"
    [n_set, n_stimuli, n_ori, _] = retrieve_dim(stim_energy)
    
    # Preallocate matrices
    mins  = np.zeros((n_set, n_ori, n_stimuli))
    means = np.zeros((n_set, n_ori, n_stimuli))
    maxs  = np.zeros((n_set, n_ori, n_stimuli))
    
    # Loop and fill matrices
    for i in range(n_set):         # e.g. Cartesian, Polar
        for h in range(n_ori):     # orientation
            for j in range(n_stimuli):  # stimulus index

                arr = stim_energy[i, j] 

                if arr.ndim == 2:
                    # Reduced case: (1, n_ori) → direct scalar
                    val = arr[0, h]
    
                    mins[i, h, j]  = val
                    means[i, h, j] = val
                    maxs[i, h, j]  = val

                    spatialrange = "entire image"
                    
                else:
                
                    img = stim_energy[i, j][representative_chIdx, h]
        
                    stats = get_stats_within_eccentricity(
                        img, pixpdeg, min_ecc_deg=4, max_ecc_deg=8
                    )
        
                    mins[i, h, j]  = stats["min"]
                    means[i, h, j] = stats["mean"]
                    maxs[i, h, j]  = stats["max"]

                    spatialrange = "within 4–8°"

    # for shared axis limits
    if metric == 'mean':
        y_all = means
    elif metric == 'max':
        y_all = maxs
    elif metric == 'min':   # note: use 'min' not 'mins' for consistency
        y_all = mins
    
    ymin = np.nanmin(y_all)
    ymax = np.nanmax(y_all)
    
    # add a small margin (optional)
    margin = 0.05 * (ymax - ymin)
    ymin -= margin
    ymax += margin
    
    # Plot results for each set (i)
    # Custom x labels for each set
    labels_set0 = ["ver", "diagR", "hor", "diagL"]
    labels_set1 = ["pin", "spirR", "ann", "spirL"]
    
    # Jitter amount (in x-axis units)
    #jitter = 0.15
    offsets = np.linspace(-0.1, 0.1, n_ori)
    
    for i in range(n_set):
        
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4))  # <-- single figure, two panels
    
        # select labels for this set
        x_labels = labels_set0 if i == 0 else labels_set1
    
        # ------------------------- SCATTER PLOT (LEFT) -------------------------
        for h in range(n_ori):
            x = np.arange(n_stimuli) + offsets[h]
    
            if metric == 'mean':
                y = means[i, h]
            elif metric == 'max':
                y = maxs[i, h]
            elif metric == 'min':   # correction
                y = mins[i, h]
    
            ax.scatter(x, y, s=40, label=ori_labels[h], alpha=0.8)
    
        # y-axis limits
        if yrange is None:
            ax.set_ylim(ymin, ymax)
        else:
            ax.set_ylim(yrange[0], yrange[1])
    
        #ax.set_title(f"Set {i} — {metric} {spatialrange}", fontsize=16)
        ax.set_title(f"{set_labels[i]} — {metric} {spatialrange}", fontsize=16)
        ax.set_xlabel("Grating type", fontsize=14)
        ax.set_ylabel(f"{metric.capitalize()} value", fontsize=14)
        ax.set_xticks(range(n_stimuli))
        ax.set_xticklabels(x_labels)
        ax.legend(title="Orientation", fontsize=8)
    
        # ------------------------- BAR PLOT (RIGHT) -------------------------
        # mean across stimuli
        if metric == 'mean':
            avg_vals = means[i].mean(axis=1)
        elif metric == 'max':
            avg_vals = maxs[i].mean(axis=1)
        elif metric == 'min':
            avg_vals = mins[i].mean(axis=1)
    
        ax2.bar(np.arange(n_ori), avg_vals, color='gray', alpha=0.8)
    
        # same y-axis limits
        if yrange is None:
            ax2.set_ylim(ymin, ymax)
        else:
            ax2.set_ylim(yrange[0], yrange[1])
    
        #ax2.set_title(f"Set {i} — avg {metric} (all gratings)", fontsize=16)
        ax2.set_title(f"{set_labels[i]} — avg {metric} (all gratings)", fontsize=16)
        ax2.set_xlabel("Orientation Channel", fontsize=14)
        ax2.set_ylabel(f"{metric.capitalize()} value", fontsize=14)
        ax2.set_xticks(np.arange(n_ori))
        ax2.set_xticklabels(ori_labels)
    
        plt.tight_layout()
        plt.show()


def canonical_normalization(stim_energy, pixpdeg, representative_chIdx=0, tuned=False):

    # note that inversion occurs with tuned=False, sigma=0.01 or 0.1, p_exp=1, q_exp=2, std_deg=3
    
    # normalization strength (when sigma is large --> less normalization; approaching 0 is strong normalization)
    sigma = 0.1
    p_exp = 1     # exponent on numerator often 1 or 2
    q_exp = 1     # The exponent q > 1 makes suppression superlinear, meaning strong signals get disproportionately more 
                  # normalized than weak ones. This matches divisive gain control models in which suppressive drive grows faster 
                  # than linear with local contrast. Excitation grows slower than suppression
    
    # std of spatial gaussian
    std_deg = 3                  # make this 1 to 5
    std_pix = std_deg * pixpdeg
    
    #std_pix = np.shape(stim_energy[1,1][1,1])[0]
    
    [n_set, n_stimuli, n_ori, n_sf] = retrieve_dim(stim_energy)

    # here the input will either already be averaged into 1 dim [i,h][1,j] or we select the dimension w/representative_chIdx
    if representative_chIdx > n_sf:
        raise ValueError(f"Error: spatial frequency channels of input are {n_sf}. Likely already reduced in a previous step.")

    # initialized output -- will be [i,h][1,j]
    norm_energy = {
    key: np.full((1,) + val.shape[1:], np.nan)
    for key, val in stim_energy.items()
    }

    suppression_field = {
    key: np.full((1,) + val.shape[1:], np.nan)
    for key, val in stim_energy.items()
    }
    
    energy = {
        key: np.full((1,) + val.shape[1:], np.nan)
        for key, val in stim_energy.items()
    }
    
    for i in range(n_set):              # e.g. Cartesian, Polar
        for j in range(n_stimuli):  # stimulus index
            for h in range(n_ori):          # orientation

                energy[i,j][0,h] = stim_energy[i, j][representative_chIdx, h]

                if tuned==False:
                    # ----------------------------------------
                    # UNTUNED NORMALIZATION
                    # ----------------------------------------
                    suppression_field[i,j][0,h] = gaussian_filter(energy[i, j][0, h] ** q_exp, sigma=std_pix)
                elif tuned==True:
                    # ----------------------------------------
                    # ORIENTATION-TUNED NORMALIZATION
                    # ----------------------------------------

                    # 1) Orientation-tuned spatial surround for orientation h
                    S_h = gaussian_filter(stim_energy[i,j][0,h] ** q_exp, sigma=std_pix)

                    # 2) Cross-orientation suppression for orientation h
                    other_oris = [o for o in range(n_ori) if o != h]
                    C_h = np.sum([stim_energy[i,j][0,o] for o in other_oris], axis=0)

                    # 3) Total suppression
                    suppression_field[i,j][0,h] = S_h + C_h

                norm_energy[i,j][0,h] = energy[i,j][0,h] ** p_exp / (sigma + suppression_field[i,j][0,h])
                
    return norm_energy


def normalization_byAnisotropy(stim_energy, pixpdeg, representative_chIdx=0):

    # normalization strength (when sigma is large --> less normalization; approaching 0 is strong normalization)
    sigma = 0.1

    [n_set, n_stimuli, n_ori, n_sf] = retrieve_dim(stim_energy)

    # here the input will either already be averaged into 1 dim [i,h][1,j] or we select the dimension w/representative_chIdx
    if representative_chIdx > n_sf:
        raise ValueError(f"Error: spatial frequency channels of input are {n_sf}. Likely already reduced in a previous step.")

    energy = {
        key: np.full((1, n_ori), np.nan)
        for key in stim_energy.keys()
    }

    std_energy = {
        key: np.full((1, 1), np.nan)
        for key in stim_energy.keys()
    }

    norm_energy = {
        key: np.full((1, n_ori), np.nan)
        for key in stim_energy.keys()
    }
    
    for i in range(n_set):              # e.g. Cartesian, Polar
        for j in range(n_stimuli):  # stimulus index
            for h in range(n_ori):          # orientation

                img = stim_energy[i, j][representative_chIdx, h]
                energy[i,j][0,h] = np.sum(img) / img.size

            # now compute the std 
            vals = energy[i,j][0, :]          # shape (4,)
            std_energy[i,j][0,0] = np.std(vals)

            norm_energy[i,j][0,:] = vals / (sigma + std_energy[i,j][0,0])

    return norm_energy


def normalization_byAnisotropy_NOA(stim_energy, pixpdeg, representative_chIdx=0):

    """
    stim_energy[(i,j)] has shape (1, n_ori, X, Y)
        (spatialFreqChannel, orientation, x, y)

    Returns norm_energy with the *same shape*.
    """
    
    p_exp = 1,
    sigma = 0.1,
    #sigma_NOA=75
    
    # std of NOA
    sigma_NOA_deg = 3                  # make this 1 to 5
    sigma_NOA = sigma_NOA_deg * pixpdeg

    
    # Extract dictionary keys
    keys = list(stim_energy.keys())
    n_set     = max(k[0] for k in keys) + 1
    n_stimuli = max(k[1] for k in keys) + 1

    # number of orientations
    sample_arr = next(iter(stim_energy.values()))
    _, n_ori, X, Y = sample_arr.shape

    # Allocate output with SAME shape
    norm_energy = {
        key: np.zeros((1, n_ori, X, Y))
        for key in stim_energy.keys()
    }

    # Loop over stimuli
    for i in range(n_set):
        for j in range(n_stimuli):

            E_full = stim_energy[i, j]        # shape (1, n_ori, X, Y)
            E_ch   = E_full[representative_chIdx]   # shape (n_ori, X, Y)

            # ---------------------------------------------------------
            # STEP 1 — orientation energy vector: average over space
            # ---------------------------------------------------------
            # Compute Eori: shape (n_ori,)
            Eori = np.mean(E_ch, axis=(1,2))

            # mean across orientations
            Ebar = np.mean(Eori)

            # ---------------------------------------------------------
            # STEP 2 — orientation variance and std
            # ---------------------------------------------------------
            diffsq = (Eori - Ebar)**2
            var_E  = np.mean(diffsq)
            std_E  = np.sqrt(var_E)

            # ---------------------------------------------------------
            # STEP 3 — scalar suppressive drive s_val
            # ---------------------------------------------------------
            #denom_NOA = (sigma_NOA**2 + std_E**2)
            #s_val = np.mean((Eori**2) / denom_NOA)

            s_val = np.mean((Eori) / std_E)

            # ---------------------------------------------------------
            # STEP 4 — output normalization,
            # applied to all pixels and orientations
            # ---------------------------------------------------------
            # E_ch: shape (n_ori, X, Y)

            R_full = E_ch ** p_exp / (s_val + sigma)   # shape (n_ori, X, Y)

            # Store output (add channel dimension back)
            norm_energy[i,j][0,:,:,:] = R_full

    return norm_energy










    

    


