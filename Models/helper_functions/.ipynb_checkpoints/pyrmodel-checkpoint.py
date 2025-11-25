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

class SteerablePyramidSF:
    def __init__(self, n_sf, n_ori, img_size, img_degrees, 
                 dtype=torch.float32, is_complex=True, 
                 downsample=False, tight_frame=False, device='cpu'):
        self.n_sf = n_sf
        self.n_ori = n_ori
        self.img_size = img_size
        self.img_degrees = img_degrees
        self.is_complex = is_complex
        self.downsample = downsample
        self.tight_frame = tight_frame
        self.dtype = dtype
        self.device = device
        self.sf_peaks = []

        # Initialize the Steerable Pyramid only once
        self.filters = SteerablePyramidFreq([self.img_size, self.img_size],
                                            height=n_sf, order=n_ori - 1,
                                            is_complex=is_complex,
                                            downsample=downsample,
                                            tight_frame=tight_frame).to(device)


    def make_empty_image(self):
        empty_image = torch.zeros((1, 1, self.img_size, self.img_size), dtype=self.dtype).to(self.device)
        return empty_image

    def reconstruct_filters(self, visualize=False):
        empty_image = self.make_empty_image()
        filter_visualization = self.filters.forward(empty_image)

        # insert a 1 in the center of each coefficient
        for k, v in self.filters.pyr_size.items():
            mid = (v[0] // 2, v[1] // 2)
            filter_visualization[k][0, 0, mid[0], mid[1]] = 1
        # ... and then reconstruct this dummy image to visualize the filter.
        filt = []
        for k in filter_visualization:
            if isinstance(k, tuple):
                filt.append(self.filters.recon_pyr(filter_visualization, [k[0]], [k[1]]))
        if visualize:
            po.imshow(filt,
                col_wrap=self.n_ori,
                title=f"len = {len(filt)}",
                vrange="indep1",
                zoom=0.5)
        return filt

    def get_filter_sf_tunings(self, to_degrees=True):
        """
        Computes the radial average of the reconstructed filters' Fourier representations.

        Returns:
            torch.Tensor: Radial average tuning curves for each spatial frequency band.
        """
        # Get the reconstructed filters from `reconstruct_filters()`
        reconList = self.reconstruct_filters(visualize=True)

        # Compute the Fourier magnitude of each reconstructed filter
        fourier_filts = [torch.fft.fftshift(torch.fft.fft2(filt).abs()) for filt in reconList]

        # Arrange filters into a structured shape
        fourier_filts = einops.pack(fourier_filts, "* c h w")[0]
        fourier_filts = einops.rearrange(fourier_filts, "(scales ori) c h w -> scales (ori c) h w",
                                         scales=self.n_sf)

        # Sum across orientations to get spatial frequency tuning
        fourier_filts = fourier_filts.sum(1, keepdim=True)

        # Compute radial average
        polar_rad = po.tools.polar_radius(fourier_filts.shape[-2:]).to(torch.int)
        nr = torch.bincount(polar_rad.ravel())

        # Normalize Fourier amplitudes
        fourier_filts = (255 * fourier_filts / fourier_filts.max()).to(torch.int)

        radial_avg = []
        for sc in range(self.n_sf):
            tbin = torch.bincount(polar_rad.ravel(), fourier_filts[sc, 0].ravel())
            radial_avg.append(tbin / nr)

        radial_avg = torch.stack(radial_avg)

        return radial_avg




    def reconstruct_optimal_stimuli(self, visualize=False):
        """
        Reconstruct the optimal stimuli for each filter in the steerable pyramid.
        Returns:
            - optimal_stimuli (list): List of optimal stimuli for each filter.
        """
        sf_preferences = self.find_sf_preference_for_each_filter()
        ori_preferences = self.find_ori_preference_for_each_filter()
        optimal_stims = {}

        for i, sf in enumerate(sf_preferences):
            for j, ori in enumerate(ori_preferences):
                optimal_stims[i, j] = self.generate_sinusoidal_grating(sf, ori, 0)

        if visualize:
            #fig, axes = plt.subplots(self.n_sf, self.n_ori, figsize=(self.n_ori, self.n_sf))
            fig, axes = plt.subplots(self.n_sf, self.n_ori, figsize=(self.n_ori, self.n_sf), dpi=400)
            for i, sf in enumerate(sf_preferences):
                for j, ori in enumerate(ori_preferences):
                    im = axes[i, j].imshow(optimal_stims[i, j], cmap='gray', vmin=-1, vmax=1)
                    sf_dg = np.round(sf/self.img_degrees,2) # below was int(sf)
                    axes[i, j].set_title(f'SF: {sf_dg}cyc/deg\nOri: {int(np.rad2deg(ori))}\u00b0', 
                                         fontsize=4)
                    axes[i, j].axis('off')
                    # Add individual colorbar to each subplot
                    #cbar = plt.colorbar(im, ax=axes[i, j], fraction=0.046, pad=0.04)
                    #cbar.ax.tick_params(labelsize=7)
            plt.tight_layout()
        return optimal_stims
    

    def find_sf_preference_for_each_filter(self, to_cpd=False):
        """
        Find the spatial frequency peak in the filter bank.
        Returns:
            - sf_peaks (list): List of spatial frequency peaks for spatial frequency filters. The spatial frequency is in the unit of cycles per image.
        """
        sf_peaks = []
        for i in range(self.n_sf):
            sf_peaks.append(self.img_size/2**(i+2))
        self.sf_peaks = sf_peaks
        if to_cpd:
            self.sf_peaks = [sf / self.img_degrees for sf in self.sf_peaks]
        return self.sf_peaks

    def find_ori_preference_for_each_filter(self):
        """
        Returns:
        - ori_peaks (list): List of orientation peaks for spatial frequency filters. The orientation is in the unit of radians.
        """
        ori_peaks = np.linspace(0, np.pi, 4, endpoint=False)[[0, 3, 2, 1]]
        return ori_peaks


    def generate_sinusoidal_grating(self, cycles, theta, phase=0):
        """
        Generate a sinusoidal grating pattern.
    
        Parameters:
            image_size (int): Size of the image (e.g., 512x512).
            cycles (float): Spatial frequency (cycles per image).
            theta (float): Orientation in radians. This means that the variation of the grating is along the direction of theta.
            phase (float): Phase offset in radians.
    
        Returns:
            np.ndarray: 2D sinusoidal grating image.
        """
        frequency = cycles / self.img_size  # Convert cycles per image to cycles per pixel
        x = np.linspace(-self.img_size//2, self.img_size//2, self.img_size)
        y = np.linspace(-self.img_size//2, self.img_size//2, self.img_size)
        X, Y = np.meshgrid(x, y)
        # Rotate coordinates
        X_theta = X * np.cos(theta) - Y * np.sin(theta)
        # Generate sinusoidal wave
        grating = np.sin(2 * np.pi * frequency * X_theta + phase)
    
        return grating

    def reshape_img(self, image):
        """
        Reshape the image to the required dimensions for the steerable pyramid.
        Assumes image is either a NumPy array or a PyTorch tensor.
        Returns:
            - torch.Tensor: Reshaped image tensor [1,1,n,n].
        """
        if isinstance(image, torch.Tensor):
            image = image.detach().cpu().numpy()  # Convert to NumPy if necessary

        # Ensure image is grayscale and has correct shape
        if len(image.shape) == 3:  # If RGB, convert to grayscale
            image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

        # Resize if needed
        if image.shape != (self.img_size, self.img_size):
            print('Resizing image to', self.img_size, 'x', self.img_size)
            image = cv2.resize(image, (self.img_size, self.img_size), interpolation=cv2.INTER_CUBIC)

        # Convert back to tensor
        image = torch.tensor(image, dtype=self.dtype).to(self.device).unsqueeze(0).unsqueeze(0)  # Shape: [1,1,H,W]
        return image
        

    def forward(self, image):
        """
        Processes an image through the steerable pyramid and returns its coefficients.
        """
        image_batch = self.reshape_img(image)
        pyr_coeffs = self.filters(image_batch)
        
        return pyr_coeffs

    def normalize_filters(self, pyr_energy, is_magnitude=False):
        """
        Normalize filter output values across sf channels using L2 normalization.
        """
        normalized_dict = {}
        if is_magnitude is True:
            i = 2
        else:
            i = 4
        for (sf, ori), item in pyr_energy.items():
            normalized_dict[sf, ori] = item / (i**sf)
        return normalized_dict

    def get_energy_of_filters(self, pyr_coeffs, to_numpy=True, reshape=True, normalize=True, to_magnitude=False):
        """
        Compute the squared magnitude (energy) of the steerable pyramid coefficients.

        Args:
            pyr_coeffs: Pyramid coefficients from the forward pass.
            to_numpy (bool): Convert results to NumPy arrays.
            reshape (bool): Reshape back to the original image size.

        Returns:
            dict: Dictionary of squared magnitude filter responses.
        """
        if to_magnitude:
            i = 1
        else:
            i = 2
        pyr_energy = {(sf, ori): (pyr_coeffs[sf, ori].abs() ** i) for sf, ori in
                    product(range(self.n_sf), range(self.n_ori))}
        if to_numpy:
            pyr_energy = {k: v.cpu().numpy() for k, v in pyr_energy.items()}
        #TODO: This needs to be fixed. The downsample is not working.
        if reshape:
            pyr_energy = {k: v.reshape((self.img_size, self.img_size)) for k, v in pyr_energy.items()}
        if normalize:
            pyr_energy = self.normalize_filters(pyr_energy, is_magnitude=to_magnitude)
        pyr_energy_numpy = np.zeros((self.n_sf, self.n_ori, self.img_size, self.img_size))
        for i in range(self.n_sf):
            for j in range(self.n_ori):
                pyr_energy_numpy[i,j] = pyr_energy[(i,j)]
        return pyr_energy_numpy

    
    def visualize_filter_magnitudes(self, pyr_mags, title=None, cbar=False, share_color_range=False, pRF_loc=None, angle_in_radians=False, max_eccentricity=4.2, scale_factor=1, save_path=None):
        """
        Visualize the filter magnitude in the steerable pyramid.
        Args:
            pyr_mags (dict): Dictionary containing magnitude responses.
        """

        fig, axes = plt.subplots(self.n_sf, self.n_ori, 
                                 figsize=(self.n_ori*scale_factor, self.n_sf*scale_factor))
        if share_color_range:
            vmin = min(image.min() for image in pyr_mags.values())  # Smallest value across all images
            vmax = max(image.max() for image in pyr_mags.values())
            norm = colors.Normalize(vmin=vmin, vmax=vmax)

        for i in range(self.n_sf):
            for j in range(self.n_ori):
                ax = axes[i, j]
                img = pyr_mags[i,j]
                ax.axis('off')
                # Handle edge cases
                if share_color_range:
                    im = ax.imshow(img, cmap='gray', norm=norm)
                else:
                    im = ax.imshow(img, cmap='gray')
                # Add individual colorbar to each subplot
                if cbar is True:
                    cbar = plt.colorbar(im, ax=ax, format="%.2e", fraction=0.046, pad=0.04)
                    cbar.ax.tick_params(labelsize=8)  # Adjust tick label size if needed
                if pRF_loc is not None:
                    if isinstance(pRF_loc[0], list) or isinstance(pRF_loc[0], tuple):
                        my_colors = ['red', 'blue', 'green', 'yellow', 'purple', 'orange']
                        for c, pRF in enumerate(pRF_loc):
                            eccentricity, angle = pRF
                            row,col = mapping.find_pRF_loc_in_image(eccentricity=eccentricity,
                                                                         angle=angle,
                                                                         angle_in_radians=angle_in_radians,
                                                                         image_size=self.img_size,
                                                                         max_eccentricity=max_eccentricity)
                            ax.scatter(col,row, color=my_colors[c], s=30, marker='o')
                            ax.text(col,row, f'{pRF}', fontsize=10, color=my_colors[c])
                    else:
                        eccentricity, angle = pRF_loc
                        row,col = mapping.find_pRF_loc_in_image(eccentricity=eccentricity,
                                                                    angle=angle,
                                                                    angle_in_radians=angle_in_radians,
                                                                image_size=self.img_size,
                                                                    max_eccentricity=max_eccentricity)
                        ax.scatter(col,row, color='red', s=30, marker='o')  # Plot pRF location

        fig.suptitle(title)
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.03, hspace=0.03)
        if save_path is not None:
            fig.savefig(save_path, bbox_inches='tight', pad_inches=0)
        plt.show()





class InputStimuli:
    def __init__(self, sfs, oris, img_size, img_degrees, referenceframe, mask, paddedBondaryPerc
                 ):
        self.sfs = sfs    # in deg
        self.oris = oris  # in deg
        self.img_size = img_size
        self.img_degrees = img_degrees
        self.reference = referenceframe
        self.mask = mask
        self.paddedBondaryPerc = paddedBondaryPerc

        # convert cycles per deg to cycles per pix
        self.freq_cpp = (self.sfs * self.img_degrees) / self.img_size
        self.ori_rad = np.deg2rad(self.oris)
        self.phase = 0

        # Create the input
        self.input = self.create_gratings_matrix()


    def create_gratings_matrix(self):

        stimuli = {}

        for i, sf in enumerate(self.freq_cpp):
            for j, ori in enumerate(self.ori_rad):
                if self.reference == 'cartesian':
                    stimuli[i, j] = self.generate_cartesian_grating(sf, ori)
                elif self.reference == 'polar':
                    stimuli[i, j] = self.generate_polar_grating(sf, ori)
                else:
                    raise ValueError("reference must be 'cartesian' or 'polar'")
    
        if self.mask:
                mask = self.makeMask(0.015*self.paddedBondaryPerc, 1*self.paddedBondaryPerc, 'soft')
                stimuli = {k: v * mask for k, v in stimuli.items()}
        
        return stimuli

    def generate_cartesian_grating(self, sf, ori):

        x = np.linspace(-self.img_size//2, self.img_size//2, self.img_size)
        y = np.linspace(-self.img_size//2, self.img_size//2, self.img_size)
        X, Y = np.meshgrid(x, y)
        # Rotate coordinates
        X_theta = X * np.cos(ori) - Y * np.sin(ori)
        # Generate sinusoidal wave
        grating = np.sin(2 * np.pi * sf * X_theta + self.phase)
        
        return grating

    def generate_polar_grating(self, sf, ori):
        x = np.linspace(-self.img_size//2, self.img_size//2, self.img_size)
        y = np.linspace(-self.img_size//2, self.img_size//2, self.img_size)
        X, Y = np.meshgrid(x, y)
        
        w_a = sf*np.pi/2 * np.cos(ori)
        w_r = sf*np.pi/2 * np.sin(ori)
        
        w_r = w_r*self.img_size; w_a = w_a*self.img_size;

        th = np.arctan2(Y, X)
        r = np.sqrt(X**2 + Y**2)
        
        w_a = np.round(w_a);
        grating = np.cos(w_r * np.log(r) + w_a * th);

        return grating

    def makeMask(self, r_inner, r_outer, edgeSetting='soft'):

        halfN = self.img_size//2
        
        # normalized coords (radius ~0.5 at edge)
        x = np.arange(-halfN, halfN) / halfN
        y = np.arange(-halfN, halfN) / halfN
        x_norm, y_norm = np.meshgrid(x, y)
    
        th_norm, r_norm = np.angle(x_norm + 1j*y_norm), np.sqrt(x_norm**2 + y_norm**2)
    
        if edgeSetting == 'hard':
            # binary annulus
            mask = ((r_norm < r_outer) & (r_norm > r_inner)).astype(float)
    
        elif edgeSetting == 'soft':
    
            soft_width = 0.05
    
            # inner ramp 0→1
            inner_ramp   = (r_norm - r_inner) / soft_width
            inner_window = 0.5 - 0.5 * np.cos(np.pi * np.clip(inner_ramp, 0, 1))
            inner_taper  = (r_norm >= r_inner) * inner_window
    
            # outer ramp 1→0
            outer_ramp   = (r_outer - r_norm) / soft_width
            outer_window = 0.5 - 0.5 * np.cos(np.pi * np.clip(outer_ramp, 0, 1))
            outer_taper  = (r_norm <= r_outer) * outer_window
    
            mask = inner_taper * outer_taper
    
        else:
            raise ValueError("edgeSetting must be 'soft' or 'hard'")
    
        return mask































        
