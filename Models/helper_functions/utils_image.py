#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import torch
import cv2
import copy
from plenoptic.simulate import SteerablePyramidFreq
import plenoptic as po
import einops
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from itertools import product
from scipy.ndimage import gaussian_filter


def stimLookupTable():

    # the order of stimuli and polar angles here should not change, otherwise verify mapping is not affected.
    
    polar_angles = [0, 45, 90, 135, 180, 225, 270, 315]

    #polar_angles = list(range(8))

    sets = {
        "Cartesian Set": ["ver", "diagR", "hor", "diagL"],
        "Polar Set": ["pin", "spirL", "ann", "spirR"]
    }
    
    analysis_types = ["Cartesian analysis", "Polar analysis"]
    
    # initialize structure
    stimulus_lookup = {
        set_name: {
            stim: {
                analysis: {} for analysis in analysis_types
            }
            for stim in stim_list
        }
        for set_name, stim_list in sets.items()
    }
    
    # fill in values
    for set_name, stim_list in sets.items():
        for stim in stim_list:
    
            # Case 1 — analysis matches set → auto-fill with stim label
            if set_name == "Cartesian Set":
                stimulus_lookup[set_name][stim]["Cartesian analysis"] = {
                    i: stim for i in polar_angles
                }
            else:
                stimulus_lookup[set_name][stim]["Polar analysis"] = {
                    i: stim for i in polar_angles
                }
    
            # Case 2 — mismatched analysis → leave empty for manual filling
            if set_name == "Cartesian Set":
                stimulus_lookup[set_name][stim]["Polar analysis"] = {
                    i: None for i in polar_angles
                }
            else:
                stimulus_lookup[set_name][stim]["Cartesian analysis"] = {
                    i: None for i in polar_angles
                }

    # add values for cases where stimulus and analysis do NOT match:
    initial_manual_vals = ["ann", "spirL", "pin", "spirR", "ann", "spirL", "pin", "spirR"]
    
    manual_vals = initial_manual_vals.copy()
    stimuli = stimulus_lookup["Cartesian Set"]

    for stim in stimuli:
        for i, val in enumerate(manual_vals):
            stimulus_lookup["Cartesian Set"][stim]["Polar analysis"][polar_angles[i]] = val
        
        # rotate list: move first element to the end
        manual_vals = manual_vals[1:] + manual_vals[:1]

    # do the same for the polar stimuli
    initial_manual_vals = ["hor", "diagR", "ver", "diagL", "hor", "diagR", "ver", "diagL"]
    manual_vals = initial_manual_vals.copy()
    stimuli = stimulus_lookup["Polar Set"]

    for stim in stimuli:
        for i, val in enumerate(manual_vals):
            stimulus_lookup["Polar Set"][stim]["Cartesian analysis"][polar_angles[i]] = val
        
        # rotate list: move first element to the end
        manual_vals = manual_vals[1:] + manual_vals[:1]


    return stimulus_lookup


def get_angle_profile(setname,
                      analysis_type,
                      target_label,
                      stimorder,
                      angles,
                      y_all):
    """
    Returns an angle-indexed vector where each value is pulled from the stimulus
    whose lookup table entry matches `target_label`.

    Parameters
    ----------

    setname : str
        e.g. "Cartesian Set" or "Polar Set".

    analysis_type : str
        e.g. "Cartesian analysis" or "Polar analysis".

    target_label : str
        e.g. "ver", "diagR", "pin", "ann" etc.

    stimorder : list of str
        Order of stimuli matching row order of y_all.
        e.g. ["ver", "diagR", "hor", "diagL"]

    angles : list of int
        Polar angle bins, e.g. [0,45,90,135,180,225,270,315]

    y_all : numpy array, shape (n_stimuli, n_angles)
        Values to extract.

    Returns
    -------
    final_vec : numpy array of shape (n_angles,)
        Final angular profile.
    """

    lookup = stimLookupTable()
    
    final_vec = []

    for ang_idx, angle in enumerate(angles):

        # matched_stim = None
        # # Search all stimuli to find which maps to the target label
        # for stim in stimorder:
        #     val = lookup[setname][stim][analysis_type][angle]
        #     if val == target_label:
        #         matched_stim = stim
        #         break
        
        matched_stim = None
        for stim in stimorder:

            # print(setname)
            # print(analysis_type)
            # print(angle)
            # print(target_label)
            
            if lookup[setname][stim][analysis_type][angle] == target_label:
                matched_stim = stim
                break

        # If nothing matches, store NaN
        if matched_stim is None:
            final_vec.append(np.nan)
            print('NO VALUE FOUND -- placing nan')
        else:
            # Retrieve row index
            stim_idx = stimorder.index(matched_stim)
            #print(stim_idx)
            final_vec.append(y_all[stim_idx, ang_idx])

    return np.array(final_vec)

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

def get_stats_within_eccentricity_by_angle(
        img, pixpdeg, min_ecc_deg=4, max_ecc_deg=8,
        angle_centers=[0,45,90,135,180,225,270,315]):
    """
    Compute min, max, mean in an eccentricity annulus AND angular bins.

    Returns:
        stats = {
            'min':  (8,) array,
            'max':  (8,) array,
            'mean': (8,) array
        }
    """
    
    # --- Convert eccentricity limits to pixels ---
    min_ecc_px = min_ecc_deg * pixpdeg
    max_ecc_px = max_ecc_deg * pixpdeg

    # --- Image dims and center ---
    ny, nx = img.shape
    cy, cx = ny/2, nx/2

    # Coord grids
    y, x = np.ogrid[:ny, :nx]

    # Radius map
    r = np.sqrt((x - cx)**2 + (y - cy)**2)

    # Eccentricity mask
    ecc_mask = (r >= min_ecc_px) & (r <= max_ecc_px)

    # --- Compute angle map (in degrees 0..360) ---
    theta = np.degrees(np.arctan2(-(y - cy), (x - cx)))  # negative y for screen coords
    theta = (theta + 360) % 360

    # Half-width of angular bin
    halfw = 22.5  # degrees

    min_vals = []
    max_vals = []
    mean_vals = []

    # --- Loop over the 8 angle bins ---
    for ac in angle_centers:
        # angular window: ac ± 22.5
        lo = (ac - halfw) % 360
        hi = (ac + halfw) % 360

        if lo < hi:
            ang_mask = (theta >= lo) & (theta < hi)
        else:
            # wrap-around case (e.g. bin centered at 0°)
            ang_mask = (theta >= lo) | (theta < hi)

        # combine radial + angular mask
        combined = ecc_mask & ang_mask

        vals = img[combined]

        if vals.size > 0:
            min_vals.append(float(vals.min()))
            max_vals.append(float(vals.max()))
            mean_vals.append(float(vals.mean()))
        else:
            min_vals.append(np.nan)
            max_vals.append(np.nan)
            mean_vals.append(np.nan)

    return {
        "min":  np.array(min_vals),
        "max":  np.array(max_vals),
        "mean": np.array(mean_vals)
    }

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


def plot_energy_per_stim(stim_energy, representative_chIdx, analysis_name, pixpdeg, min_ecc_deg=4, max_ecc_deg=8, yrange=None):

    polar_angles = [0, 45, 90, 135, 180, 225, 270, 315]
    n_pa = len(polar_angles)
    
    #ori_labels = ["90°", "45°", "0°", "135°"]

    stimorder_set0 = ["ver", "diagR", "hor", "diagL"]
    stimorder_set1 = ["pin", "spirR", "ann", "spirL"]
    
    set_labels = ["Cartesian Set", "Polar Set"]

    analyses = ["Cartesian analysis", "Polar analysis"]
    n_analyses = len(analyses)
    
    metric="mean"
    [n_set, n_stimuli, n_ori, _] = retrieve_dim(stim_energy)
    
    # Preallocate matrices
    mins  = np.zeros((n_set, n_stimuli, n_pa))
    means = np.zeros((n_set, n_stimuli, n_pa))
    maxs  = np.zeros((n_set, n_stimuli, n_pa))
    
    # Loop and fill matrices
    for i in range(n_set):                  # e.g. Cartesian, Polar Stimulus set
        for j in range(n_stimuli):          # stimulus index

            # sum the energy/magnitude across orientation channels
            img = np.squeeze(np.sum(stim_energy[i,j], axis=1))
            
            # for h in range(n_pa):           # polar angles (8)
            #     for k in range(n_analyses): # Cartesian or Polar analysis
        
            stats = get_stats_within_eccentricity_by_angle(
                img, pixpdeg, min_ecc_deg=4, max_ecc_deg=8
            )
        
            mins[i, j, :]  = stats["min"]
            means[i, j, :] = stats["mean"]
            maxs[i, j, :]  = stats["max"]

            spatialrange = "within 4–8°"

    # for shared axis limits
    if metric == 'mean':
        y_all = means
    elif metric == 'max':
        y_all = maxs
    elif metric == 'min':   # note: use 'min' not 'mins' for consistency
        y_all = mins

    if analysis_name == "Cartesian analysis": ############ WHICH ANALYSIS
        targetlist = stimorder_set0
        analysis_name = analyses[0] # do cartesian only for now
        colorpro = [0.498, 0.749, 0.4824]
        colorcon = [0.6863, 0.5529, 0.7647]
        colorneut = [0.5, 0.5, 0.5]         
        
    elif analysis_name == "Polar analysis":
        targetlist = stimorder_set1
        analysis_name = analyses[1] # do cartesian only for now

        colorpro = [0.5725, 0.7725, 0.8706]
        colorcon = [0.7922, 0, 0.1255]
        colorneut = [0.5, 0.5, 0.5] 
    
    for setname in set_labels:                ############ WHICH STIMULUS SET
        
        print(setname)
        set_idx = set_labels.index(setname)     # 0 or 1
        currData = np.squeeze(y_all[set_idx, :, :])

        if setname == "Cartesian Set":
            stimorder = stimorder_set0
        elif setname == "Polar Set":
            stimorder = stimorder_set1


        angles_deg = np.array([0,45,90,135,180,225,270,315])
        angles_rad = np.deg2rad(angles_deg)
        
        fig = plt.figure(figsize=(12,6))
        
        # ---- LEFT: POLAR PLOT ----
        ax_polar = fig.add_subplot(1, 2, 1, polar=True)

        
        for idx, target_label in enumerate(targetlist):   # e.g., ["ver","diagR","hor","diagL"]
            
            profile = get_angle_profile(
                setname=setname,
                analysis_type=analysis_name,
                target_label=target_label,
                stimorder=stimorder,
                angles=angles_deg,
                y_all=currData
            )
        
            # close circular loop
            profile_closed = np.append(profile, profile[0])
            angles_closed  = np.append(angles_rad, angles_rad[0])
        
            # same color logic
            if analysis_name == "Cartesian analysis": 
                if idx in [0, 2]:
                    color = colorpro
                else:
                    color = colorcon
            elif analysis_name == "Polar analysis": 
                if idx in [0]:
                    color = colorpro
                elif idx in [2]:
                    color = colorcon
                else:
                    color = colorneut
        
            ax_polar.plot(angles_closed, profile_closed, marker='o', linewidth=2,
                          label=target_label, color=color)
    
        ax_polar.legend(loc='upper right')
        ax_polar.set_ylim(yrange[0], yrange[1])
        ax_polar.set_title("Local direction energy", pad=20)
        ax_polar.set_theta_zero_location("E")
        ax_polar.set_theta_direction(1)
        
        # ---- RIGHT: BAR PLOT ----
        ax_bar = fig.add_subplot(1, 2, 2)
        
        avg_values = []
        bar_colors = []
        
        for idx, target_label in enumerate(targetlist):
        
            profile = get_angle_profile(
                setname=setname,
                analysis_type=analysis_name,
                target_label=target_label,
                stimorder=stimorder,
                angles=angles_deg,
                y_all=currData
            )
        
            avg_values.append(np.nanmean(profile))

            if analysis_name == "Cartesian analysis": 
                if idx in [0, 2]:
                    bar_colors.append(colorpro)
                else:
                    bar_colors.append(colorcon)
            elif analysis_name == "Polar analysis": 
                if idx in [0]:
                    bar_colors.append(colorpro)
                elif idx in [2]:
                    bar_colors.append(colorcon)
                else:
                    bar_colors.append(colorneut)
        
        ax_bar.bar(targetlist, avg_values, color=bar_colors)
        ax_bar.set_ylabel("Average energy")
        ax_bar.set_title("Mean across polar angles")
        
        if yrange is not None:
            ax_bar.set_ylim(yrange[0], yrange[1])
        
        
        plt.tight_layout()
        plt.show()


def canonical_normalization(stim_energy, pixpdeg, representative_chIdx=0, p_exp = 1, q_exp = 1, tuned=False):

    # note that inversion occurs with tuned=False, sigma=0.01 or 0.1, p_exp=1, q_exp=2, std_deg=3
    
    # normalization strength (when sigma is large --> less normalization; approaching 0 is strong normalization)
    sigma = 0.1
    
    # p_exp: exponent on numerator often 1 or 2
    # q_exp: The exponent q > 1 makes suppression superlinear, meaning strong signals get disproportionately more 
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

            # this will be used to compute the overall suppressive drive accounting for all orientations
            pooled = np.sum(stim_energy[i,j][0,:] ** q_exp, axis=0)     # pool across orientations first
            #np.shape(stim_energy)
            Z = gaussian_filter(pooled, sigma=std_pix)             # convolve pooled energy across all orientation channels
            
            for h in range(n_ori):          # orientation

                energy[i,j][0,h] = stim_energy[i, j][representative_chIdx, h]

                if tuned==False:
                    # ----------------------------------------
                    # UNTUNED NORMALIZATION
                    # ----------------------------------------

                    # use this if want to compute suppression per orientation channel
                    #suppression_field[i,j][0,h] = gaussian_filter(energy[i, j][0, h] ** q_exp, sigma=std_pix)
                    #print(np.shape(suppression_field[i,j][0,h]))
                    #print(type(suppression_field[i,j][0,h]))

                    # use this if want to use the same suppressive drive for all orientation channels
                    suppression_field[i,j][0,h] = Z
                    
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





def div_normalization(stim_energy, pixpdeg, p_exp = 1, q_exp = 1, tuned=False):
    
    # normalization strength (when sigma is large --> less normalization; approaching 0 is strong normalization)
    sigma = 0.01
    
    # p_exp: exponent on numerator often 1 or 2
    # q_exp: The exponent q > 1 makes suppression superlinear, meaning strong signals get disproportionately more 
                  # normalized than weak ones. This matches divisive gain control models in which suppressive drive grows faster 
                  # than linear with local contrast. Excitation grows slower than suppression
    
    # std of spatial gaussian
    std_deg = 3                  # make this 1 to 5
    std_pix = std_deg * pixpdeg
    
    [n_set, n_stimuli, n_ori, n_sf] = retrieve_dim(stim_energy)

    # initialized output -- will be [i,h][1,j]
    norm_energy = {
        key: np.full((n_sf,) + val.shape[1:], np.nan)
        for key, val in stim_energy.items()
    }

    Z = {
        key: np.full((n_sf,) + val.shape[1:], np.nan)
        for key, val in stim_energy.items()
    }
    
    for i in range(n_set):              # e.g. Cartesian, Polar
        for j in range(n_stimuli):  # stimulus index

            if tuned==False:

                # ----------------------------------------
                # UNTUNED NORMALIZATION
                # ----------------------------------------
                
                # apply 4D gaussian filter across 4D matrix (all SFs, ORIs, X, Y)
                Z[i,j] = gaussian_filter(
                    stim_energy[i, j] ** q_exp,
                    sigma=(std_pix, std_pix, std_pix, std_pix)
                )

            elif tuned==True:

                # ----------------------------------------
                # ORIENTATION-TUNED NORMALIZATION
                # ----------------------------------------
            
                for h in range(n_ori):          # orientation

                    # 1) Orientation-tuned spatial surround for orientation h ((6, 1, 896, 896))
                    S_h = gaussian_filter(
                        stim_energy[i,j][:,h:h+1] ** q_exp, 
                        sigma=(std_pix, 1, std_pix, std_pix)
                    )

                    # 2) Cross-orientation suppression for orientation h
                    other_oris = [o for o in range(n_ori) if o != h]

                    C_h = np.sum(
                        [stim_energy[i,j][:, o:o+1] for o in other_oris],  # each (6, 1, 896, 896)
                        axis=0
                    ) 

                    # 3) Total suppression
                    Z[i,j][:,h:h+1] = S_h + C_h

            norm_energy[i,j] = stim_energy[i, j] ** p_exp / (sigma + Z[i,j])
            
    return norm_energy






def normalization_byAnisotropy(stim_energy, pixpdeg):

    # normalization strength (when sigma is large --> less normalization; approaching 0 is strong normalization)
    sigma = 0.01
    
    # take mean across SF channels
    stim_energy_SFave = {
        key: np.mean(val, axis=0, keepdims=True)
        for key, val in stim_energy.items()
    }

    [n_set, n_stimuli, n_ori, n_sf] = retrieve_dim(stim_energy_SFave)

    energy = {
        key: np.full((n_sf, n_ori), np.nan)
        for key in stim_energy_SFave.keys()
    }

    std_energy = {
        key: np.full((1, 1), np.nan)
        for key in stim_energy_SFave.keys()
    }

    norm_energy = {
        key: np.full((n_sf, n_ori), np.nan)
        for key in stim_energy_SFave.keys()
    }
    
    for i in range(n_set):              # e.g. Cartesian, Polar
        for j in range(n_stimuli):      # stimulus index
            for s in range(n_sf):      # sf
                for h in range(n_ori):      # orientation

                    img = stim_energy_SFave[i, j][s, h]
                    energy[i,j][s,h] = np.sum(img) / img.size

                # now compute the std 
                vals = energy[i,j][s, :]          # shape (4,)
                
                std_energy[i,j][s,0] = np.std(vals)

                norm_energy[i,j][s,:] = vals / (sigma + std_energy[i,j][s,0])

    return norm_energy


def normalization_byAnisotropy_NOA(stim_energy, pixpdeg):

    # STD in denominator is NON-SPATIAL

    """
    stim_energy[(i,j)] has shape (1, n_ori, X, Y)
        (spatialFreqChannel, orientation, x, y)

    Returns norm_energy with the *same shape*.
    """
    
    p_exp = 1,
    sigma = 0.01,
    #sigma_NOA=75
    
    # std of NOA
    #sigma_NOA_deg = 3                  # make this 1 to 5
    #sigma_NOA = sigma_NOA_deg * pixpdeg

    # number of orientations
    sample_arr = next(iter(stim_energy.values()))
    _, _, X, Y = sample_arr.shape

    # take mean across SF channels
    stim_energy_SFave = {
        key: np.mean(val, axis=0, keepdims=True)
        for key, val in stim_energy.items()
    }

    [n_set, n_stimuli, n_ori, n_sf] = retrieve_dim(stim_energy_SFave)
    
    # Allocate output with SAME shape
    norm_energy = {
        key: np.zeros((n_sf, n_ori, X, Y))  
        for key in stim_energy_SFave.keys()
    }

    # Loop over stimuli
    for i in range(n_set):
        for j in range(n_stimuli):

            E_full = stim_energy_SFave[i, j]        # shape (1, n_ori, X, Y)

            # orientation energy vector: average over space
            # Compute Eori: shape (n_ori,)
            Eori = np.mean(E_full, axis=(2,3))

            # mean across orientations
            Ebar = np.mean(Eori)

            # orientation variance and std
            diffsq = (Eori - Ebar)**2
            var_E  = np.mean(diffsq)
            std_E  = np.sqrt(var_E)

            # scalar suppressive drive s_val
            s_val = std_E
            
            # E_ch: shape (n_ori, X, Y)
            R_full = E_full ** p_exp / (s_val + sigma)   # shape (n_ori, X, Y)

            # Store output (add channel dimension back)
            norm_energy[i,j] = R_full

    return norm_energy


def normalization_byStimHomogeneity(stim_energy, pixpdeg):

    # Normalization based on homogeneity of center-surround stimulus match
    # The measure of homogeneity is contrast-independent, 
    
    # Homogeneity is defined as the cosine similarity between locally pooled center and surround feature-energy vectors in the 24-D 
    # SF×orientation space. This homogeneity value is then used to weight the exponent from 1-2 in the suppressive drive.

    """
    stim_energy[(i,j)] has shape (1, n_ori, X, Y)
        (spatialFreqChannel, orientation, x, y)

    Returns norm_energy with the *same shape*.
    """

    sigma = 0.01
    p_exp = 1
    sigma_center_deg = 3
    sigma_surround_deg = 10
    
    # Convert degrees → pixels
    sigma_c = sigma_center_deg * pixpdeg
    sigma_s = sigma_surround_deg * pixpdeg
    
    # average across SF?
    # stim_energy_SFave = {
    #     key: np.mean(val, axis=0, keepdims=True)
    #     for key, val in stim_energy.items()
    # }

    # for a given stimulus,
    # go through each ori channel-- if the local and surround match ori give higher weight of suppression q (q 1-->2)
    # sum contrast of this match

    [n_set, n_stimuli, n_ori, n_sf] = retrieve_dim(stim_energy)
    # number of orientations
    sample_arr = next(iter(stim_energy.values()))
    _, _, X, Y = sample_arr.shape

    # Prepare output dict with same keys (q_exp == homogeneity, from 1-2, 2 being identical center-surround distributions)
    q_exp = {
        key: np.full(stim_energy[key].shape[2:], np.nan)
        for key in stim_energy.keys()
    }

    Z = {
        key: np.full(stim_energy[key].shape[2:], np.nan)
        for key in stim_energy.keys()
    }

    norm_energy = {
        key: np.zeros((n_sf, n_ori, X, Y))   # key: np.zeros((n_sf, n_ori, X, Y))
        for key in stim_energy.keys()
    }

    for i in range(n_set):
        for j in range(n_stimuli):

            E_full = stim_energy[i, j]        # (6, 4, 896, 896)

            ######## CALCULATE HOMOGENEITY ACROSS SPACE) ##########
            # --- Center and surround via convolution --- 
            C = gaussian_filter(E_full, sigma=(0, 0, sigma_c, sigma_c)) # weighted average at center neighborhood at x,y
            S = gaussian_filter(E_full, sigma=(0, 0, sigma_s, sigma_s)) # weighted average for larger pool (surround) at x,y

            # To normalize distributions:
            # --- Compute sum across SF and ORI dimensions ---
            C_sum = np.sum(C, axis=(0, 1), keepdims=True) + 1e-9  # total center energy across ALL 24 channels (per x,y)
            S_sum = np.sum(S, axis=(0, 1), keepdims=True) + 1e-9  # total center energy across ALL 24 channels (per x,y)

            # --- Divide the weighted average (local / surround) by their respective sums 
            # to enable the magnitude comparisons between local and surround (converts raw energies into channel distributions that sum to 1)
            # This makes the homogeneity measure contrast-independent
            P_c = C / C_sum
            P_s = S / S_sum

            # --- Flatten channel dimension: (24, H, W) from (6, 4, H, W) ---
            Pc_flat = P_c.reshape(-1, *P_c.shape[2:])
            Ps_flat = P_s.reshape(-1, *P_s.shape[2:])

            ########### COSINE SIMILARITY PER PIXEL ############
            # compute the dot product between the 24-D vectors at each pixel 
            # ("match scores" for center vs. surround on a channel by channel basis)
            dot = np.sum(Pc_flat * Ps_flat, axis=0)      # (H,W)

            # compute overall strength of center / surround to normalize the dot
            norm_c = np.sqrt(np.sum(Pc_flat**2, axis=0))
            norm_s = np.sqrt(np.sum(Ps_flat**2, axis=0))

            # this division makes homogeneity NOT depend on selectivity (uniform channels can still be homogenous)
            # cos(angle) = (A · B) / (||A|| · ||B||)
            H_map = dot / (norm_c * norm_s + 1e-9)
            #######################################################
            
            # Store (0 ---> very different channel distribution; 1 --> very similar channel distribution between center/surround)
            q_exp[i, j] = H_map + 1 # added plus one because this will be used for the q exponent, which should range from 1-2 not 0-1

            # apply 4D gaussian filter across 4D matrix (all SFs, ORIs, X, Y) -- exponent scaled based on surround similarity
            Z[i,j] = gaussian_filter(
                stim_energy[i, j] ** q_exp[i,j],
                sigma=(sigma_c, sigma_c, sigma_c, sigma_c)
                )

            norm_energy[i,j] = stim_energy[i, j] ** p_exp / (sigma + Z[i,j])
            
    return norm_energy


def condense_energy(energy_dict, dim=0, operation="mean"):
    """
    energy_dict[n_set, n_stim] -> array (D0, D1, H, W)
    Returns same structure, but reduced along `dim` with size=1.
    """

    out = {}

    for key, arr in energy_dict.items():
        if operation == "mean":
            reduced = np.mean(arr, axis=dim, keepdims=True)
        elif operation == "sum":
            reduced = np.sum(arr, axis=dim, keepdims=True)
        else:
            raise ValueError("operation must be 'mean' or 'sum'")

        out[key] = reduced

    return out


def orientation_mask(theta, ori_angles, h, w=22.5):
    centers = ori_angles[h]
    mask = np.zeros_like(theta, dtype=bool)

    for ang in centers:
        d = np.abs((theta - ang + 180) % 360 - 180)  # circular distance
        mask |= (d <= w)

    return mask

def orientation_weight(theta, ori_angles, h):
    """
    theta: (H, W) array of angles in degrees, in [0, 360)
    ori_angles: dict like {0: [0, 180], 1: [45, 225], ...}
    h: orientation channel index

    Returns: (H, W) weight array in [0, 1]
             1 at preferred directions (e.g. 90, 270),
             0 at orthogonal (±90° away, e.g. 0, 180).
    """
    centers = ori_angles[h]              # e.g. [90, 270] for vertical
    # start with large distances
    d_min = np.full_like(theta, np.inf, dtype=float)

    for ang in centers:
        # circular distance in degrees, in [0, 180]
        d = np.abs((theta - ang + 180) % 360 - 180)
        d_min = np.minimum(d_min, d)

    # map distance 0 → 1, distance 90° → 0, clip beyond 90 to 0
    w = 1.0 - np.clip(d_min / 90.0, 0.0, 1.0)

    return w
    

def imposeImbalance(stim_energy, cardinal=1, radial=1):

    w = 22.5

    [n_set, n_stimuli, n_ori, n_sf] = retrieve_dim(stim_energy)
    sample_arr = next(iter(stim_energy.values()))
    _, _, W, H = sample_arr.shape
    
    yc, xc = (H - 1) / 2, (W - 1) / 2

    y, x = np.meshgrid(np.arange(H), np.arange(W), indexing='ij')
    theta = (np.degrees(np.arctan2(y - yc, x - xc)) + 360) % 360

    stim_energy_imb = copy.deepcopy(stim_energy)
    
    
    for i in range(n_set):             # e.g. Cartesian, Polar
        for j in range(n_stimuli):     # orientation

            ###### CARDINAL OVERREPRESENTATION
            ## do regardless of whether cartesian or polar b/c cardinal channels are always in cartesian categories (0 & 2)
            # vertical channel
            if cardinal!=1:
                stim_energy_imb[i, j][:, 0] *= cardinal
                # horizontal channel
                stim_energy_imb[i, j][:, 2] *= cardinal

            ##### RADIAL OVERREPRESENTATION
            if radial!=1:
                # for polar, need to define radial based on polar angle and channel
                ori_angles = {
                    0: [90, 270], # vertical channel (radial is 90/270)
                    1: [135,315], # up-right channel (radial is 45/225)
                    2: [0, 180],  # horizontal channel (radial is 0/180)
                    3: [45, 225] # up-left channel (radial is 135/315)
                }
    
                #ori_masks = {h: orientation_mask(theta, ori_angles, h, w) for h in range(4)}
                ori_masks = {h: orientation_weight(theta, ori_angles, h) for h in range(4)}
                
                # for h in range(4):
                #     stim_energy_imb[i, j][:, h][:, ori_masks[h]] *= radial
    
                for h in range(4):
                    w_h = ori_masks[h]          # (896, 896), values in [0, 1]
                    N = radial                     # desired max gain
                
                    # scale from 1 (no change) up to N at preferred directions
                    gain = 1.0 + (N - 1.0) * w_h   # (896, 896), in [1, 2]
                
                    stim_energy_imb[i, j][:, h] *= gain  # broadcast over SF dimension
                

    return stim_energy_imb


    

    


