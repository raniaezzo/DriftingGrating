% tinker vignetting

clear; 

addpath(genpath(fullfile(pwd, 'SimplifiedSteerablePyramid')))
stimtype = 'polar';
%stimtype = 'cartesian';

numorientations = 4;
numscales = 6;

ppd = 512; % simulate as entire width of image as 1 deg, then adjust frequency etc.

%frequencies = ppd ./ 2.^(1:numscales+2);
setIdx = 4;                                         % set idx 4 to match exp stimulus

% let's manualy define frequencies so they match actual stimulus
experimenStimDiam = 24;
base = experimenStimDiam ./ 2.^(1:numscales+2);     % exponentially space the freqs
k = 1 / base(setIdx);                                    % scale so 4th element = 1 (stimulus freq)
vals = k * base;                                    % these are the values relative to stimulus
frequencies = experimenStimDiam * vals;                   % these account for pixel scaling

if strcmp(stimtype, 'polar')
    % 4 because for radial freq to match at midpoint this is diameter/4
    frequencies = frequencies; %frequencies*pi*0.5;
end

whichscale = setIdx;

% these will be scaled (here 1 deg is diam; in exp 12 deg is diam)
% so 1 cycle per deg must be * by 12
stimfrequency = frequencies(whichscale); % cycles per degree

orientations = pi*(0:numorientations-1)/numorientations; % in radians
[x, y] = meshgrid((-ppd/2:ppd/2-1)/(ppd/2)); % make this ppd/2 instead of ppd to double padding
[th, r] = cart2pol(x,y);

% this should represent the actual stimulus edges (not analysis region)
%mask = double(r<.4 & r >.1);

% simulated DIAMETER is 1 deg so divide answer by 2 to scale
% inner cut out = 0.34 deg (std of gauss) / 12.2 max ~= 0.03 deg / 2
% outer cut out = 12.2 deg / 12.2 max ~= 1 deg / 2 
mask = double(r<0.5 & r >.015);



for stim_ori = 1:numorientations
    
    switch stimtype
        case 'cartesian'
            kx = stimfrequency*sin(orientations(stim_ori));
            ky = stimfrequency*cos(orientations(stim_ori));

            stim = cos(2*pi*(kx*x+ky*y));
            IM(:,:,stim_ori) = stim .* mask;
        case 'polar'
             % tried (stimfrequency/(1/pi))/sqrt(2)
            [w_a, w_r] = pol2cart(orientations(stim_ori), stimfrequency); %/sqrt(2)); % was 10      
            w_r = w_r; w_a = w_a;   % *4 to match the half radius (6/24)?
            w_a = round(w_a);
            im = cos(w_r * log(r) + w_a * th) .* mask;
            im(isnan(im))=0;
            IM(:,:,stim_ori) = im;
    end

    coeff = buildSCFpyr(IM(:,:,stim_ori),numscales+2, numorientations);
    coeff = coeff{whichscale};
    sz = size(coeff{1}); 
    coeff = reshape(cell2mat(coeff), sz(1), sz(2), []);

    E(:,:,stim_ori) = sum(abs(coeff),3);

    

end

figure(1); clf; 
montage(IM, DisplayRange=[-1 1])
title('input images')
figure(2); clf
montage(E);
title('energy per input image (across SF and Ori channels')
figure(3); clf
montage(E-mean(E,3), DisplayRange=4*[-1 1]);
title('mean-centered energy per input image')

%%
% RE: to get the min/max vals
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

global_max = 83.5; global_min = 0.0000001; % just set a constant

% steerable pyramid gives different size outputs per freq channel by design 
selectFreqChOutput = zeros(128, 128, numorientations, numorientations); % (x, y, orientation, input)

for s = 1:numorientations %input
    figure(4+s-1),
    coeff = buildSCFpyr(IM(:,:,s),numscales+2, numorientations);
    tiledlayout(numscales,numorientations)
    for ii = 1:numscales %SF
        for jj=1:numorientations %OR
            nexttile();
            vals = abs(coeff{ii+1}{jj}).^2;
            %scaled_vals = (vals - global_min) ./ (global_max - global_min);
            scaled_vals = (vals - global_min) ./ (global_max - global_min);

            imshow(scaled_vals)
            %imshow((abs(coeff{ii+1}{jj}).^2));
            
            if ii==4 % used to be 3 but i changed the SF range
                scaled_all(:, :, jj, s) = scaled_vals;
            end
        end
    end
    sgtitle(sprintf('ori channel %i', s))
end


%% Mask out anything not in the analysis region
%     % maybe this should be done AFTER normalization?
% 
% % Assume scaled_all is 128x128x4x4
% sz = size(scaled_all);
% H = sz(1);
% W = sz(2);
% 
% % Define center (assuming centered grid)
% cx = (W+1)/2;
% cy = (H+1)/2;
% 
% % Compute radius map
% [x, y] = meshgrid(1:W, 1:H);
% r = sqrt((x - cx).^2 + (y - cy).^2);
% 
% % Define inner and outer radius thresholds
% r_inner = H * (1/3);
% r_outer = H * (2/3);
% 
% % Build mask: keep only pixels between inner and outer radius
% mask = (r >= r_inner) & (r <= r_outer);
% 
% % Apply to all 4x4 slices
% masked_all = scaled_all .* mask;   % MATLAB will broadcast mask to all slices


%%

% scaled_all is 128 x 128 x 4 x 4 (x, y, orientation, input)

% first account for imbalanced population of the ver/hor channels
% imbalanced neuronal sample as a gain param
weighted_energy(:, :, 1, :) = scaled_all(:, :, 1, :) * 2;
weighted_energy(:, :, 3, :) = scaled_all(:, :, 3, :) * 2;

% normalize so min is 0 and max is 1 across all units and channels
% Get min and max across first 3 dimensions for each element in 4th dim
%normalized_wenergy = rescale(scaled_all, 0, 1);
normalized_wenergy = rescale(weighted_energy, 0, 1);
%normalized_wenergy = weighted_energy;

% sum across orientation channels (no divisive normalization)
popResponse = popresponse_minmax(normalized_wenergy, 'no normalization');

% DENOMINATOR

% 1. first just apply good old normalization (same on all channels)
% this will not produce oblique > cardinal orientation
% just should just balance out the responses despite unequal units

% iteratively create gaussian size based on the weight matrix
gaussian_filter = create_gaussFilter(weighted_energy, .1);

% convolve the normalized_wenergy with this gaussian
% this convolves each 2D energy image for (:,:, n, m)
conv_wenergy = imfilter(normalized_wenergy, gaussian_filter, ...
    'replicate', 'same');

% sum across channels
conv_wenergy_acrossCh = squeeze(sum(conv_wenergy, 3));

% plot energy field across channels, for each input image
figure
subplot(2,2,1)
imshow(popResponse(:,:,1))
colorbar;
hold on; draw_center_circles(popResponse)
subplot(2,2,2)
imshow(popResponse(:,:,2))
colorbar;
hold on; draw_center_circles(popResponse)
subplot(2,2,3)
imshow(popResponse(:,:,3))
colorbar;
hold on; draw_center_circles(popResponse)
subplot(2,2,4)
imshow(popResponse(:,:,4))
colorbar;
hold on; draw_center_circles(popResponse)
sgtitle('energy E across ori channels')

% plot suppression field across channels, for each input image
figure
subplot(2,2,1)
imshow(conv_wenergy_acrossCh(:,:,1))
hold on; draw_center_circles(conv_wenergy_acrossCh)
subplot(2,2,2)
imshow(conv_wenergy_acrossCh(:,:,2))
hold on; draw_center_circles(conv_wenergy_acrossCh)
subplot(2,2,3)
imshow(conv_wenergy_acrossCh(:,:,3))
hold on; draw_center_circles(conv_wenergy_acrossCh)
subplot(2,2,4)
imshow(conv_wenergy_acrossCh(:,:,4))
hold on; draw_center_circles(conv_wenergy_acrossCh)
sgtitle('energy normalizer D across ori channels')

figure
const = 0.01; % higher number means less suppression
% lower temporal frequency channel:
subplot(2,2,1)
imshow(popResponse(:,:,1)./(conv_wenergy_acrossCh(:,:,1)+const), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(popResponse)
subplot(2,2,2)
imshow(popResponse(:,:,2)./(conv_wenergy_acrossCh(:,:,2)+const), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(popResponse)
subplot(2,2,3)
imshow(popResponse(:,:,3)./(conv_wenergy_acrossCh(:,:,3)+const), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(popResponse)
subplot(2,2,4)
imshow(popResponse(:,:,4)./(conv_wenergy_acrossCh(:,:,4)+const), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(popResponse)
sgtitle('Normalized contrast energy across ori channels (more suppression)')

figure
const = 0.1; % higher number means less suppression
% high temporal frequency channel:
subplot(2,2,1)
imshow(popResponse(:,:,1)./(conv_wenergy_acrossCh(:,:,1)+const), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(popResponse)
subplot(2,2,2)
imshow(popResponse(:,:,2)./(conv_wenergy_acrossCh(:,:,2)+const), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(popResponse)
subplot(2,2,3)
imshow(popResponse(:,:,3)./(conv_wenergy_acrossCh(:,:,3)+const), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(popResponse)
subplot(2,2,4)
imshow(popResponse(:,:,4)./(conv_wenergy_acrossCh(:,:,4)+const), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(popResponse)
sgtitle('Normalized contrast energy across ori channels (less suppression)')

%% 2. for tuned normalization

% basically pixel unit
gaussian_filter_null = create_gaussFilter(weighted_energy, .01);

% convolve the normalized_wenergy with this gaussian
% this convolves each 2D energy image for (:,:, n, m)

tunedNorm_resp = nan(size(normalized_wenergy));
const = 0.01; %0.01; % higher number means less suppression (or 0.1)

conv_wenergy_tuned = zeros(size(normalized_wenergy), 'like', normalized_wenergy);
for k = 1:size(normalized_wenergy,3)
    for m = 1:size(normalized_wenergy,4)
        if k == m
            conv_wenergy_tuned(:,:,k,m) = imfilter(normalized_wenergy(:,:,k,m), gaussian_filter, 'replicate', 'same');
        else
            conv_wenergy_tuned(:,:,k,m) = imfilter(normalized_wenergy(:,:,k,m), gaussian_filter_null, 'replicate', 'same');
        end
        % do the channel-specific normalization
        tunedNorm_resp(:,:,k,m) = normalized_wenergy(:,:,k,m) ./ (conv_wenergy_tuned(:,:,k,m) + const);
    end
end

tunedNorm_allCh = sum(tunedNorm_resp,3) ./ size(tunedNorm_resp,3);

figure
subplot(2,2,1)
imshow(tunedNorm_allCh(:,:,1), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(tunedNorm_allCh)
subplot(2,2,2)
imshow(tunedNorm_allCh(:,:,2), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(tunedNorm_allCh)
subplot(2,2,3)
imshow(tunedNorm_allCh(:,:,3), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(tunedNorm_allCh)
subplot(2,2,4)
imshow(tunedNorm_allCh(:,:,4), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(tunedNorm_allCh)
sgtitle('Tuned Normalized contrast energy across ori channels')

%% 2. for tuned normalization version 2

% basically pixel unit
gaussian_filter_null = create_gaussFilter(weighted_energy, .01);

% convolve the normalized_wenergy with this gaussian
% this convolves each 2D energy image for (:,:, n, m)

tunedNorm_resp = nan(size(normalized_wenergy));
const = 0.01; %0.01; % higher number means less suppression (or 0.1)


conv_wenergy_tuned = zeros(size(normalized_wenergy), 'like', normalized_wenergy);

for m = 1:size(normalized_wenergy,4) % inputs
    for k = 1:size(normalized_wenergy,3) % channels
        
        currentCh = imfilter(normalized_wenergy(:,:,k,m), gaussian_filter, 'replicate', 'same');
        
        otherChs = imfilter(normalized_wenergy(:,:,k,m), gaussian_filter_null, 'replicate', 'same');
        
    end

    %conv_wenergy_tuned(:,:,k,m) = 

    % do the channel-specific normalization
    tunedNorm_resp(:,:,k,m) = normalized_wenergy(:,:,k,m) ./ (conv_wenergy_tuned(:,:,k,m) + const);
end

tunedNorm_allCh = sum(tunedNorm_resp,3) ./ size(tunedNorm_resp,3);

figure
subplot(2,2,1)
imshow(tunedNorm_allCh(:,:,1), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(tunedNorm_allCh)
subplot(2,2,2)
imshow(tunedNorm_allCh(:,:,2), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(tunedNorm_allCh)
subplot(2,2,3)
imshow(tunedNorm_allCh(:,:,3), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(tunedNorm_allCh)
subplot(2,2,4)
imshow(tunedNorm_allCh(:,:,4), [], 'DisplayRange', [0 1])
hold on; draw_center_circles(tunedNorm_allCh)
sgtitle('Tuned Normalized contrast energy across ori channels')


%%
% create a normalization weight matrix to operate on weighted_energy 

% --- maybe I can just use a convolution for each weighted energy channel:
% GAUSSIAN SCALED BY ITS ENERGY (Energy 0 --> .01 , and Energy 1 --> .1)
    % if the energy channel == image orientation, apply large gaus (.1)
    % if the energy channel ~= image orientation, apply pixel gaus (.01)

% for each channel, do this..? yes! if there is a match, .1
    % not sure this will work for pinwheel, but let's try



%% functions

function popResponse_norm = popresponse_minmax(scaled_all, titlename)
    popResponse = squeeze(sum(scaled_all, 3));

    % then divide by the # of channels to keep min/max
    popResponse_norm = popResponse./size(scaled_all,3);
    %popResponse_norm = popResponse;

    % % then normalize: this channel did not have a max of 1, due to lower SF 
    % % around the edge -- make this the case b/c I am only analyzing this:
    % % Normalize the entire array to [0, 1] <-- assumes perfect match and
    % % perfect mismatch are present in this array (min and max)
    % min_val = min(popResponse(:));
    % max_val = max(popResponse(:));
    % popResponse_norm = (popResponse - min_val) / (max_val - min_val);
    
    % figure
    % % before suppression (higher count)
    % subplot(2,2,1)
    % imshow(popResponse_norm(:,:,1))
    % subplot(2,2,2)
    % imshow(popResponse_norm(:,:,2))
    % subplot(2,2,3)
    % imshow(popResponse_norm(:,:,3))
    % subplot(2,2,4)
    % imshow(popResponse_norm(:,:,4))
    % sgtitle(titlename)
end

function gaussian_filter = create_gaussFilter(channelImages, sigma2imageFraction)
    % create the standard gaussian filter 1/10th of the image size (24 deg/10).
    [H, W] = size(channelImages, 1, 2); % Get image size
    sigma_h = H * sigma2imageFraction;
    sigma_w = W * sigma2imageFraction;
    
    % Define filter size to be 6*sigma (common practice), rounded to nearest odd
    filter_size_h = 2 * floor(3 * sigma_h) + 1;
    filter_size_w = 2 * floor(3 * sigma_w) + 1;
    
    % Create Gaussian filter (sum of elements = 1)
    gaussian_filter = fspecial('gaussian', [filter_size_h, filter_size_w], mean([sigma_h, sigma_w]));
end 

function draw_center_circles(inputIm)
% draw_center_circles Draws red circles at 1/3 and 2/3 of image size.
%
%   [cx, cy, r1, r2] = draw_center_circles(tunedNorm_allCh)
%
%   Inputs:
%       tunedNorm_allCh : 2D or 3D array (image or stack)
%                         Circles are drawn based on size of the first two dims.
%
%   Outputs:
%       cx, cy : coordinates of image center
%       r1, r2 : radii (1/3 and 2/3 of image height)

    % Image size
    [H, W, ~] = size(inputIm);

    % Center coordinates
    cx = W / 2;
    cy = H / 2;

    % Radii (fractions of height)
    r1 = round(H * ((1/3)/2));
    r2 = round(H * ((2/3)/2));

    centers = [cx cy; cx, cy];
    radii = [r1; r2];

    % Draw red circles
    viscircles(centers, radii, 'Color', 'r', 'LineWidth', 1);

    hold off
end