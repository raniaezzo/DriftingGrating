clear;

addpath(genpath(fullfile(pwd, 'SimplifiedSteerablePyramid')))

stimtype = 'cartesian'; % cartesian or polar

numorientations = 4;   % number of orientation channels in pyramid, 
                            % also used for number of stimuli passed in
numscales = 6;         % number of bandpass scales
whichscale = 4;        % scale index you care about for later extraction
                        % coeff_all{1}   % highpass residual
                        % coeff_all{2}   % band k = 1
                        % coeff_all{3}   % band k = 2
                        % coeff_all{4}   % band k = 3  ← actual matched band??
                        % coeff_all{5}   % band k = 4  ← ??
                        % coeff_all{6}   % band k = 5
                        % coeff_all{7}   % band k = 6
                        % coeff_all{8}   % lowpass residual

%% --- IMAGE GRID SETUP (now 1024 x 1024) to include padding ---

N = 768;    % image size
halfN = N/2;

% Pixel coordinate grids centered at 0
% xx_pix, yy_pix are in pixel units (e.g. -512 ... +511)
[x_pix, y_pix] = meshgrid((-halfN:(halfN-1)), (-halfN:(halfN-1)));

% to make this in degrees (less straight forward w/ padding
%[x_pix, y_pix] = meshgrid((-halfN:(halfN-1))/halfN);

% Make the outer radius half of max (rest will be padding)
r_inner = 0.015; r_outer = 0.5;
mask = makeMask(halfN, r_inner, r_outer, 'soft'); % soft or hard edges


%% --- GRATINGS WITH INTEGER NUMBER OF CYCLES ACROSS THE IMAGE ---

% Ideally an integer number of total cycles across the screen diameter.

% Pick something reasonable, not too close to Nyquist.
nCyclesAcross = 48; % double it since I shrunk the stimulus to half for padding 
                    % this is how many cycles across the unmasked image

base_freq_cycPerPix = nCyclesAcross / N; % cycles/pixel magnitude

% orientations to test (0, 45, 90, 135 deg)
stim_oris = pi * (0:numorientations-1) / numorientations;  % [0, pi/4, pi/2, 3pi/4]

IM = zeros(N, N, numorientations);  % store each oriented grating (masked)

for oi = 1:numorientations

    thetaV = stim_oris(oi);

    switch stimtype
        case 'cartesian'
            % Frequency vector in cycles/pixel
            fx = base_freq_cycPerPix * cos(thetaV);
            fy = base_freq_cycPerPix * sin(thetaV);
        
            % Cartesian grating: cos(2π (fx*x + fy*y))
            stim = cos( 2*pi * (fx * x_pix + fy * y_pix) );
        
            % Apply circular mask
            stim_masked = stim .* mask;
        
            IM(:,:,oi) = stim_masked;
        case 'polar'
            [w_a, w_r] = pol2cart(thetaV, base_freq_cycPerPix);
            w_r = w_r*N; w_a = w_a*N;

            [th, r] = cart2pol(x_pix,y_pix);

            w_a = round(w_a);
            im = cos(w_r * log(r) + w_a * th) .* mask;
            im(isnan(im))=0;
            IM(:,:,oi) = im;
    end
end

% quick check that each input image is matched in power
rms_val = @(img) sqrt(mean(img(:).^2));
for oi = 1:4
    raw_rms(oi) = rms_val(IM(:,:,oi));
end
raw_rms

%% --- PYRAMID DECOMPOSITION AND ENERGY MAPS ---

E = zeros(N, N, numorientations); % energy maps for each input image

for oi = 1:numorientations

    coeff_all = buildSCFpyr(IM(:,:,oi), numscales+2, numorientations);
    % coeff_all is {highpassResidual}{band1}{band2}...{bandN}{lowpassResidual}
    % bandpass scales are 2 (first and last?).
    % Your "whichscale" should refer to one of those band cells.

    this_scale_cell = coeff_all{whichscale};  % <- scale you care about
    % this_scale_cell is 1 x numorientations, each complex subband

    % Stack orientation channels at this spatial scale
    sz_sub = size(this_scale_cell{1});
    coeff_mat = reshape(cell2mat(this_scale_cell), sz_sub(1), sz_sub(2), []);

    % Sum energy across orientation channels at that spatial frequency band
    E(:,:,oi) = imresize( sum(abs(coeff_mat),3), [N N], "nearest" );
    % note: coeff_mat is smaller than N x N for deeper scales.
    % we upscale back to N×N for visualization/comparison.
end

figure(1); clf;
montage(IM, 'DisplayRange', [-1 1]);
title(sprintf('Input gratings (%ix%i, masked, integer cycles)', N, N));

figure(2); clf;
montage(E, 'DisplayRange', []);
title('Energy per input image at chosen scale');



% what is the image size of each level?

fprintf('scale-  ii-  subbandSize-downsampleFactor-cycPerPix-cycAcrossIm\n')
for ii = 2:(numscales+1)  % skip highpass at coeff_all{1}
    % Take one orientation subband (they all have same size)
    sz = size(coeff_all{ii}{1});
    N_band = sz(1);  % assume square
    downsample_factor = N / N_band;

    % Estimate center frequency based on downsampling ratio
    % Nyquist frequency = 0.5 cyc/pix
    f_cycPerPix = 0.5 / downsample_factor;

    fprintf('%5d | %2d | %4d x %-4d | %8.1f× | %10.5f | %10.1f\n', ...
        ii-1, ii, N_band, N_band, downsample_factor, f_cycPerPix, f_cycPerPix*N);
end

% check the scale (SF level) to use in the analysis

for ii = 2:numscales+1
    band_s = coeff_all{ii};
    subcat = cat(3, band_s{:});
    band_energy(ii) = sum(abs(subcat(:)).^2);  % use squared magnitude (power)
end
[~, best_idx] = max(band_energy);
fprintf('Best ii = %d (coeff_all{%d})\n', best_idx, best_idx);

%% --- TOTAL PYRAMID ENERGY ACROSS ALL SCALES / ORIENTATIONS ---

pyr_energy = zeros(1,numorientations);

for oi = 1:numorientations
    % this returns complex-values output coeff_all{s}{o}
    % where {s} indexes the spatial freq band, and {o} the orientation
    coeff_all = buildSCFpyr(IM(:,:,oi), numscales+2, numorientations);

    total_energy = 0; % across SFs
    for s = 2:(numscales+1)             % bandpass scales (skip HF)
        band_s = coeff_all{s};           % 1 x numorientations
        subcat = cat(3, band_s{:});      % stack orientations

        % how strongly image drives orientation channels at that SF
        mag = abs(subcat).^2;               % complex magnitude
                                            % square to compute as 'power'

        total_energy = total_energy + sum(mag(:)); % collapse across pixels x,y, and ori chs
    end
    pyr_energy(oi) = total_energy;
end

disp('Total steerable energy per orientation (all bandpass scales pooled):');
disp(pyr_energy);

% each input image total energy
ratio_horiz_vs_vert = pyr_energy(1) / pyr_energy(3); % assuming ori 1=0 rad, ori 3=pi/2
fprintf('Horizontal / Vertical energy ratio = %.4f\n', ratio_horiz_vs_vert);

%%

global_min = Inf;
global_max = -Inf;

for s = 1:numorientations
    coeff = buildSCFpyr(IM(:,:,s),numscales+2, numorientations);

    for ii = 1:numscales
        for jj = 1:numorientations
            vals = abs(coeff{ii+1}{jj}).^2;
            current_min = min(vals(:));
            current_max = max(vals(:));
            if current_min < global_min
                global_min = current_min;
            end
            if current_max > global_max
                global_max = current_max;
                ind = ii;
            end
        end
    end
end
%


%% KEEP IN MIND: stimulus order was changed to match channel order

global_max = 208.6248; global_min = 0; %.0000001; % just set a constant

% steerable pyramid gives different size outputs per freq channel by design 
scaled_all = zeros(192, 192, numorientations, numorientations); % (x, y, orientation, input)

if strcmp(stimtype, 'cartesian')
    stim_ori_labels = {'ver', 'upR', 'hor', 'upL'};
elseif strcmp(stimtype, 'polar')
    stim_ori_labels = {'pin', 'spCC', 'ann', 'spC'};
end

channel_labels = {'ver', 'upR', 'hor', 'upL'};

for s = 1:numorientations %input
    figure(4+s-1),
    coeff = buildSCFpyr(IM(:,:,s),numscales+2, numorientations);
    tiledlayout(numscales,numorientations)
    for ii = 2:numscales+1 %SF
            for jj=1:numorientations %OR
                nexttile();
                vals = abs(coeff{ii}{jj}).^2;
                scaled_vals = (vals - global_min) ./ (global_max - global_min);
                %scaled_vals = vals;
    
                % this is just to visualize everything normalized by the
                % same min/max: the channel of a freq dominates
                % interesting that cross comparison is possible across SF
                % channels given how the model is described
                %imshow(scaled_vals) % keep this for now

                % imagesc(scaled_vals)
                clim = [global_min global_max];
                imagesc(vals, clim);
                axis square
                set(gca,'xtick',[]); set(gca,'ytick',[])
                colormap gray
                
                if ii==whichscale % when the scales matches the freq
                    % note that I save out vals, which is not scaled
                    % and I verified that these values are reasonable
                    % given the commented block directly above
                    scaled_all(:, :, jj, s) = vals; % scaled_vals; %
                end
            end
    end
    sgtitle(sprintf('input image %s (ch: ver, upR., hor, upL)', stim_ori_labels{s}))
end

%%

% at one frequency:
% scaled_all is 128 x 128 x 4 x 4 (x, y, orientation, input)

% Sum across x and y dimensions
valsSummary = squeeze(sum(sum(scaled_all, 1), 2));  % result is 4x4 (orientation × input)


% Alternatively, visualize as 2D matrix (orientations × inputs)
figure;
clims = [0 1300000];
%clims = [0 1165000];
imagesc(valsSummary, clims);
hold on;
n = size(valsSummary, 1);

% Draw grid lines between cells
for i = 0.5 : 1 : n+0.5
    % Horizontal lines
    plot([0.5, n+0.5], [i, i], 'k-', 'LineWidth', 1);
    % Vertical lines
    plot([i, i], [0.5, n+0.5], 'k-', 'LineWidth', 1);
end
hold off;

colorbar;
xlabel('Input orientation');
ylabel('Orientation channel');
title('energy per orientation and input');
xticks(1:4);
xticklabels(stim_ori_labels);
yticks(1:4);
yticklabels(channel_labels);


%%

% scaled_all is 128 x 128 x 4 x 4 (x, y, orientation, input)

% first account for imbalanced population of the ver/hor channels
% imbalanced neuronal sample as a gain param
weighted_energy(:, :, 1, :) = scaled_all(:, :, 1, :) .* 2;
weighted_energy(:, :, 2, :) = scaled_all(:, :, 2, :) .* 1;
weighted_energy(:, :, 3, :) = scaled_all(:, :, 3, :) .* 2;
weighted_energy(:, :, 4, :) = scaled_all(:, :, 4, :) .* 1;

% Sum across x and y dimensions
vals = squeeze(sum(sum(weighted_energy, 1), 2));  % result is 4x4 (orientation × input)

figure;
tiledlayout(2,2); % 4 subplots arranged in 2x2 grid

for ch = 1:4
    nexttile;
    bar(vals(ch,:), 'FaceColor', [0.3 0.6 0.9]);
    %ylim([0 1]);
    title(['Channel: ' channel_labels{ch}]);
    xticks(1:4);
    xticklabels(stim_ori_labels);
    xlabel('Input orientation');
    ylabel('energy');
    xlim([0 5])
    ylim([0 1300000*2])
end

%%
% Alternatively, visualize as 2D matrix (orientations × inputs)
figure;
clims = [0 1300000*2];
%clims = [0 1165000*2];
imagesc(vals, clims);
hold on;
n = size(vals, 1);

% Draw grid lines between cells
for i = 0.5 : 1 : n+0.5
    % Horizontal lines
    plot([0.5, n+0.5], [i, i], 'k-', 'LineWidth', 1);
    % Vertical lines
    plot([i, i], [0.5, n+0.5], 'k-', 'LineWidth', 1);
end
hold off;

colorbar;
xlabel('Input orientation');
ylabel('Orientation channel');
title('energy per orientation and input (2x more cardinal neurons)');
xticks(1:4);
xticklabels(stim_ori_labels);
yticks(1:4);
yticklabels(channel_labels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature-tuned normalization (original method - testing 1)  [NO MATCH]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = weighted_energy;  % [Nx, Ny, nOri, nStim]

% ---------------- Hyperparameters ----------------
sigma_space_wide  = 75;    % px, for matched orientation/stimulus
sigma_space_narrow = 0.5;  % px, ~pixel-local for non-matched
p          = 1;            % exponent on numerator
sigma0     = 0.1;          % semi-saturation constant

tunedNorm_resp = zeros(size(E), 'like', E);

for m = 1:size(E,4)        % stimulus index
    for k = 1:size(E,3)    % orientation channel index

        this_map = E(:,:,k,m);

        if k == m
            % matched orientation: broad surround suppression
            pooled = imfilter(this_map, sigma_space_wide, 'replicate', 'same');
        else
            % non-matched: essentially pixelwise/no surround
            pooled = imfilter(this_map, sigma_space_narrow, 'replicate', 'same');
            % (or just pooled = this_map; if you want literally none)
        end

        tunedNorm_resp(:,:,k,m) = this_map ./ (pooled + sigma0);
    end
end

axesMax = [2, 40000];
titlestr='feature-tuned normalization';
plot_4stimOutput(tunedNorm_resp, titlestr, stim_ori_labels, axesMax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature-tuned normalization (original method - testing 2)  [NO MATCH]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = weighted_energy;  % [Nx, Ny, nOri, nStim]

% ---------------- Hyperparameters ----------------
sigma_space_wide  = 75;    % px, for matched orientation/stimulus
sigma_space_narrow = 0.5;  % px, ~pixel-local for non-matched
p          = 1;            % exponent on numerator
sigma0     = 0.1;          % semi-saturation constant

% ---------------- Sizes ----------------
[Nx, Ny, nOri, nStim] = size(E);

% ---------------- Preallocate ----------------
S_simple   = zeros(Nx, Ny, nOri, nStim);  % surround / spatial pool
C_simple   = zeros(Nx, Ny, nOri, nStim);  % cross-orientation pool (no blur)
Z_simple   = zeros(Nx, Ny, nOri, nStim);  % denominator
R_simple   = zeros(Nx, Ny, nOri, nStim);  % final normalized output

% =================================================
% Main loop
% =================================================
for s = 1:nStim

    % ---- (A) Spatial surround, orientation-specific ----
    for o = 1:nOri
        this_ori_map = E(:,:,o,s);

        if o == s
            % matched: broad spatial pooling
            S_simple(:,:,o,s) = imgaussfilt(this_ori_map, sigma_space_wide);
        else
            % non-matched: essentially no spatial pooling
            % sigma ~ 0 makes Gaussian almost identity, but imgaussfilt
            % will complain at exactly 0, so use a tiny sigma
            S_simple(:,:,o,s) = imgaussfilt(this_ori_map, sigma_space_narrow);
        end
    end

    % ---- (B) Cross-orientation: sum other orientation channels locally ----
    for o0 = 1:nOri
        others = setdiff(1:nOri, o0);
        C_simple(:,:,o0,s) = sum(E(:,:,others,s), 3);
    end

    % ---- (C) Denominator pool ----
    Z_simple(:,:,:,s) = S_simple(:,:,:,s) + C_simple(:,:,:,s);

    % ---- (D) Divisive normalization ----
    R_simple(:,:,:,s) = ( E(:,:,:,s) .^ p ) ./ ( Z_simple(:,:,:,s) + sigma0 );
end

axesMax = [6.5, 25000];
titlestr='feature-tuned normalization';
plot_4stimOutput(R_simple, titlestr, stim_ori_labels, axesMax)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature-tuned normalization (Fang et al. style?) [NO MATCH]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = weighted_energy;  % [Nx, Ny, nOri, nStim]

[Nx, Ny, nOri, nStim] = size(E);

% Hyperparameters
                    % define in pixels/768 (but in deg is about 5)
sigma_space = 75;   % spatial surround std in *pixels* 
p            = 1;  % exponent on numerator (often 1 or 2)
sigma0       = 0.1; % semi-saturation constant (>0)

% Preallocate
S = zeros(Nx, Ny, nOri, nStim);   % same-orientation surround term
C = zeros(Nx, Ny, nOri, nStim);   % cross-orientation local term
Z_tuned   = zeros(Nx, Ny, nOri, nStim);
R_tuned   = zeros(Nx, Ny, nOri, nStim);

for s = 1:nStim

    % ---- (A) Orientation-tuned surround: blur each orientation map ----
    for o = 1:nOri
        this_ori_map = E(:,:,o,s);                           % Nx x Ny

        % weighted sum
        S(:,:,o,s)   = imgaussfilt(this_ori_map, sigma_space); % Gaussian blur in space

        % % sum within circle
        % r = 75;  % radius in pixels (you choose this)
        % % Make a disk kernel
        % [x, y] = meshgrid(-r:r, -r:r);
        % diskMask = (x.^2 + y.^2) <= r.^2;  % 1s inside the circle, 0s outside
        % % Convolve your map with that disk
        % S(:,:,o,s)   =  conv2(this_ori_map, double(diskMask), 'same');
    end

    % ---- (B) Cross-orientation suppression: sum other orientations at same pixel ----
    for o0 = 1:nOri
        others = setdiff(1:nOri, o0);
        C(:,:,o0,s) = sum(E(:,:,others,s), 3);  % no blur here
    end

    % ---- Combine to get the tuned normalization pool ----
    Z_tuned(:,:,:,s) = S(:,:,:,s) + C(:,:,:,s);

    % ---- Divisive normalization with exponent p ----
    R_tuned(:,:,:,s) = ( E(:,:,:,s) .^ p ) ./ ( Z_tuned(:,:,:,s) + sigma0 );

end

axesMax = [6.5, 25000];
titlestr='feature-tuned normalization';
plot_4stimOutput(R_tuned, titlestr,stim_ori_labels, axesMax)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature-tuned normalization (suppression superlinear in energy) [MATCH!]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = weighted_energy;  % [Nx, Ny, nOri, nStim]

% if q = p, then suppression and excitation grow at the same rate

q = 2;          % superlinear pooling exponent (>1)
                % The exponent q > 1 makes suppression superlinear, 
                % meaning strong signals get disproportionately more 
                % normalized than weak ones. This matches divisive gain 
                % control models in which suppressive drive grows faster 
                % than linear with local contrast.
p = 1;          % excitation grows slower than suppression
sigma_space = 75;
sigma0 = 0.1;

for s = 1:nStim
    for o = 1:nOri
        this_ori_map = E(:,:,o,s);  % Nx x Ny

        % surround from SAME orientation, but with superlinear drive
        S(:,:,o,s) = imgaussfilt( this_ori_map .^ q, sigma_space );
    end

    % cross-orientation
    for o0 = 1:nOri
        others = setdiff(1:nOri, o0);
        C(:,:,o0,s) = sum(E(:,:,others,s), 3);
    end

    Z_tuned(:,:,:,s) = S(:,:,:,s) + C(:,:,:,s);
    
    % clim = [0 max(Z_tuned(:,:,1,s), [], 'all')];
    % imagesc(Z_tuned(:,:,1,s), clim);
    % axis square
    % set(gca,'xtick',[]); set(gca,'ytick',[])
    % colormap gray

    R_tuned(:,:,:,s) = ( E(:,:,:,s) .^ p ) ./ ( Z_tuned(:,:,:,s) + sigma0 );
end

axesMax = [0.3, 2000];
titlestr='feature-tuned normalization (superlinear)';
plot_4stimOutput(R_tuned, titlestr, stim_ori_labels, axesMax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature-tuned normalization (by orientation anisotropy) [NO MATCH]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = weighted_energy;  % [Nx, Ny, nOri, nStim]

% Hyperparameters
p        = 1;      % exponent on numerator
sigma0   = 0.1;    % semi-saturation const in final normalization
sigma_NOA = 75;     % <-- this is your NOA sigma

[Nx, Ny, nOri, nStim] = size(E);

R_NOA = zeros(Nx, Ny, nOri, nStim);

for s = 1:nStim

    % ---- 1. Compute oriented energy per channel (spatial average) ----
    % Eori will be nOri x 1
    Eori = squeeze( mean( mean( E(:,:,:,s), 1 ), 2 ) );  % avg over x,y

    % mean across orientations
    Ebar = mean(Eori);

    % ---- 2. Compute orientation spread (std across orientations) ----
    % use population std (1/nOri) for model stability
    diffsq = (Eori - Ebar).^2;                 % nOri x 1
    var_E  = mean(diffsq);                     % scalar
    std_E  = sqrt(var_E);                      % scalar

    % ---- 3. Build scalar suppressive drive s for this stimulus ----
    % numerator: square Eori to match units of var_E
    % denominator: sigma_NOA^2 + std_E^2 ( = sigma_NOA^2 + var_E )
    denom_NOA = (sigma_NOA.^2 + std_E.^2);     % scalar

    s_val = mean( (Eori.^2) ./ denom_NOA );    % scalar, 1/nOri * sum_theta ...

    % ---- 4. Apply divisive normalization using this global s_val ----
    % same s_val across all pixels and orientations for this stimulus
    for o = 1:nOri
        R_NOA(:,:,o,s) = ( E(:,:,o,s) .^ p ) ./ ( s_val + sigma0 );
    end

end

axesMax = [3500, 13000000];
titlestr='feature-tuned normalization (global anisotropy model)';
plot_4stimOutput(R_NOA, titlestr, stim_ori_labels, axesMax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model with imbalanced population math [MATCH!]
% No spatial suppression component.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = scaled_all;

% ---- Setup ----
[Nx, Ny, nOri, nStim] = size(E);

ori_labels = rad2deg(stim_oris);   % must match channel order in E
sigma_norm = 30;                    % same as in your population code
sigma_const = 0.5;                  % same role as sigma_const in pop model

% adjust channel weights
p_balanced = [0.25  0.25  0.25  0.25 ];
p_biased   = [0.375 0.125 0.375 0.125];

% ---- Build orientation interaction matrix ----
circDiff180 = @(a,b) abs(mod(a - b + 90, 180) - 90);

W_ori = zeros(nOri, nOri);
for o1 = 1:nOri
    for o2 = 1:nOri
        d = circDiff180(ori_labels(o1), ori_labels(o2));
        W_ori(o1,o2) = exp(-(d.^2)/(2*sigma_norm^2));
    end
end

% Scale columns by population proportions
W_bal_ori = W_ori .* reshape(p_balanced, [1 nOri]);
W_bia_ori = W_ori .* reshape(p_biased,   [1 nOri]);

% ---- Allocate outputs ----
R_bal_img   = zeros(Nx, Ny, nOri, nStim);
R_bia_img   = zeros(Nx, Ny, nOri, nStim);

% ---- Per stimulus, per pixel, per orientation ----
for s = 1:nStim
    for x = 1:Nx
        for y = 1:Ny
            
            v = squeeze(E(x,y,:,s));  % nOri x 1 energy at this pixel

            % balanced pool:
            % den_bal(o1) = sigma_const + sum_o2 W_bal_ori(o1,o2)*v(o2)
            den_bal = sigma_const + (W_bal_ori * v);    % nOri x 1
            den_bia = sigma_const + (W_bia_ori * v);    % nOri x 1

            % normalized response:
            R_bal_img(x,y,:,s) = v ./ den_bal;
            R_bia_img(x,y,:,s) = v ./ den_bia;
        end
    end
end

axesMax = [7, 65000];
titlestr='balanced population implementation';
plot_4stimOutput(R_bal_img, titlestr, stim_ori_labels,axesMax)
titlestr='imbalanced population implementation';
plot_4stimOutput(R_bia_img, titlestr, stim_ori_labels,axesMax)


%%

% functions

function mask = makeMask(halfN, r_inner, r_outer, edgeSetting)

    % Normalized coordinates for mask (radius ~0.5 at edge like before)
    % match previous behavior where r<0.5 was "outer edge"
    [x_norm, y_norm] = meshgrid((-halfN:(halfN-1))/halfN, (-halfN:(halfN-1))/halfN);
    [th_norm, r_norm] = cart2pol(x_norm, y_norm);

    if strcmp(edgeSetting, 'hard')
        % Circular annulus mask (binary), scaled similarly to your old r<0.5 & r>0.015
        mask = double(r_norm < 0.5 & r_norm > 0.015);  % same proportions, now on 1024x1024
    
    elseif strcmp(edgeSetting, 'soft')
    
        soft_width = 0.05;  % fractional width of the cosine transition
        
        % Initialize
        mask_soft = zeros(size(r_norm));
        
        % Inner cosine ramp: goes from 0 → 1 between r_inner and r_inner+soft_width
        inner_ramp = (r_norm - r_inner) / soft_width;
        inner_window = 0.5 - 0.5*cos(pi*min(max(inner_ramp,0),1));  % raised cosine
        inner_taper = double(r_norm >= r_inner) .* inner_window + double(r_norm < r_inner)*0;
        
        % Outer cosine ramp: goes from 1 → 0 between r_outer-soft_width and r_outer
        outer_ramp = (r_outer - r_norm) / soft_width;
        outer_window = 0.5 - 0.5*cos(pi*min(max(outer_ramp,0),1));
        outer_taper = double(r_norm <= r_outer) .* outer_window + double(r_norm > r_outer)*0;
        
        % Combine inner and outer tapers to make a soft annulus
        mask = inner_taper .* outer_taper;
    
    end
end

% calculating output of the 4 inputs:

function plot_4stimOutput(R_tuned, titlestr, stim_ori_labels, axesMax)


    % simple sum
    R_sum = squeeze( sum(R_tuned, 3) );  
    % size: [Nx, Ny, nStim]
    
    % root of sum of squres (energy?)
    R_rss = squeeze( sqrt( sum(R_tuned.^2, 3) ) );
    % size: [Nx, Ny, nStim]
    
    
    % Compute global min and max across all channels
    clims = [0, axesMax(1)]; %max(R_tuned, [], 'all')];
    max(R_tuned, [], 'all')

    figure;
    montage(R_sum, 'Size', [2 2], 'DisplayRange', clims);
    colormap gray;
    axis image off;
    title(titlestr);
    set(gcf, 'Color', 'w');
    
    
    % Alternatively, visualize as 2D matrix (orientations × inputs)
    % Sum across x and y dimensions
    eSummary = squeeze(sum(sum(R_sum, 1), 2));  % result is 4x4 (orientation × input)
    eSummary = [eSummary(1:2)';eSummary(3:4)'];
    max(eSummary, [], 'all')

    % figure;
    % clims = [0 axesMax(2)]; %max(eSummary, [], 'all')];
    % imagesc(eSummary, clims);
    % hold on;
    % n = size(eSummary, 1);
    % title(titlestr);
    % colorbar;

    figure;
    bar([eSummary(1,:), eSummary(2,:)], 'k');  % flatten matrix column-wise into a single vector
    title(titlestr);
    ylabel('Response (a.u.)');
    ylabel('Energy across channels');
    clims = [0, axesMax(2)]; %max(eSummary, [], 'all')];
    ylim(clims);
    xlabel('Input orientation');
    xticks(1:4);
    xticklabels(stim_ori_labels);
    set(gcf, 'Color', 'w');
end