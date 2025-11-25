% implement steerable pyramid for cartesian and polar gratings
% then implement normalization

clear; 

addpath(genpath(fullfile(pwd, 'SimplifiedSteerablePyramid')))
%stimtype = 'polar';
stimtype = 'cartesian';

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
            %mask = ones(512);
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

    % this returns the numscales+2 (highpass, level1, 2, 3, 4, 5, 6,
    % lowpass)
    coeff = buildSCFpyr(IM(:,:,stim_ori),numscales+2, numorientations);

    % ignore the extra high pass by skipping
    coeff = coeff{whichscale}; %+1};
    sz = size(coeff{1}); 

    % input image N that include real/imaginary component
    % for each orientation channel response
    coeff = reshape(cell2mat(coeff), sz(1), sz(2), []);

    % energy across channels for a given input image
    E(:,:,stim_ori) = sum(abs(coeff),3);

    

end

figure(1); clf; 
montage(IM, DisplayRange=[-1 1])
title('input images')

% these will be ordered based on the filters in SCF
figure(2); clf
montage(E, DisplayRange=20*[-1 1]);
title('energy per input image (across SF and Ori channels')
%figure(3); clf
%montage(E-mean(E,3), DisplayRange=4*[-1 1]);
%title('mean-centered energy per input image')

%%

% RE: to get the min/max vals
% global_min = Inf;
% global_max = -Inf;
% 
% for s = 1:numorientations
%     coeff = buildSCFpyr(IM(:,:,s),numscales+2, numorientations);
% 
%     for ii = 1:numscales
%         for jj = 1:numorientations
%             vals = abs(coeff{ii+1}{jj}).^2;
%             current_min = min(vals(:));
%             current_max = max(vals(:));
%             if current_min < global_min
%                 global_min = current_min;
%             end
%             if current_max > global_max
%                 global_max = current_max;
%                 ind = ii;
%             end
%         end
%     end
% end
% %

global_max = 83.5; global_min = 0.0000001; % just set a constant

% steerable pyramid gives different size outputs per freq channel by design 
scaled_all = zeros(128, 128, numorientations, numorientations); % (x, y, orientation, input)

stim_ori_labels = {'hor', 'upR', 'ver', 'upL'};
channel_labels = {'ver', 'upR', 'hor', 'upL'};

for s = 1:numorientations %input
    figure(4+s-1),
    coeff = buildSCFpyr(IM(:,:,s),numscales+2, numorientations);
    tiledlayout(numscales,numorientations)
    for ii = 1:numscales %SF
        if ii==1
            continue
        else
            for jj=1:numorientations %OR
                nexttile();
                vals = abs(coeff{ii}{jj}).^2; % used to be ii+1 but now the image is mostly background
                scaled_vals = (vals - global_min) ./ (global_max - global_min);
                %scaled_vals = vals;
    
                imshow(scaled_vals)
                %imshow((abs(coeff{ii+1}{jj}).^2));
                
                if ii==setIdx % when the scales matches the freq
                    scaled_all(:, :, jj, s) = scaled_vals;
                end
            end
        end
    end
    sgtitle(sprintf('input image %s (ch: ver, upR., hor, upL)', stim_ori_labels{s}))
end

%%

% at one frequency:
% scaled_all is 128 x 128 x 4 x 4 (x, y, orientation, input)

% Sum across x and y dimensions
vals = squeeze(sum(sum(scaled_all, 1), 2));  % result is 4x4 (orientation × input)

% Optional: verify size
size(vals)   % should print [4 4]

% Flatten into 16-element vector for simple plotting
vals_flat = vals(:);

% Plot as bar graph
figure;
bar(vals_flat);
xlabel('Orientation–Input combination');
ylabel('Summed response');
title('Summed energy across orientations and inputs');

% Alternatively, visualize as 2D matrix (orientations × inputs)
figure;
imagesc(vals);
colorbar;
xlabel('Input index');
ylabel('Orientation index');
title('Summed energy per orientation and input');
xticks(1:4);
xticklabels(stim_ori_labels);
yticks(1:4);
yticklabels(channel_labels);