% I don't think it's coded for da yet.. check this
project = 'dg';

numColors = 128;
% Define RGB values for yellow, orange, and red
yellow = [1, 1, 0];
orange = [1, 0.65, 0];
red = [1, 0, 0];

% Create an array to store the RGB values
rgbArray = zeros(numColors, 3);

% Interpolate between yellow and orange
for i = 1:round(numColors/2)
    t = (i-1) / (round(numColors/2)-1); % Normalized position between 0 and 1
    rgbArray(i, :) = yellow * (1-t) + orange * t;
end

% Interpolate between orange and red
for i = round(numColors/2)+1:numColors
    t = (i-round(numColors/2)-1) / (numColors-round(numColors/2)-1); % Normalized position between 0 and 1
    rgbArray(i, :) = orange * (1-t) + red * t;
end

% for circles
% for circle
radius = 1;
center_x = 0;
center_y = 0;
theta = linspace(0, 2*pi, 100); % Define the range of theta from 0 to 2*pi
% Calculate the x and y coordinates of the circle
a = radius * cos(theta) + center_x;
b = radius * sin(theta) + center_y;
a2 = 12 * cos(theta) + center_x;
b2 = 12 * sin(theta) + center_y;


% contrasts
% contrastnames = {'cardMsep', 'oblMsep', 'allmValls', ...
%     'allsVblank', 'allmVblank', 'cardmVblank', 'oblmVblank', ...
% 8 -     'm0_v_s90', 
% 9 - 'm90_v_s0',
% 10 - 'm180_v_s90',
% 11 - 'm270_v_s0', ...
% 12 -    'm45_v_s135',
% 13 - 'm135_v_s45',
% 14 - 'm225_v_s135',
% 15 - 'm315_v_s45', 
% 16 - 'cardsVblank', 
% 17 - 'oblsVblank', ...
% 18-     'm0_v_b',
% 19- 'm180_v_b',  
% 20 - 'm90_v_b', 
% 21 - 'm270_v_b', 
% 22 - 'm45_v_b', 
% 23 - 'm225_v_b', 
% 24 - 'm135_v_b', ...
% 25 - 'm315_v_b', 
% 26 - 's0_v_b',  
% 27 - 's90_v_b', 
% 28 - 's45_v_b', 
% 29 - 's135_v_b'};

% polar angle
%polarAngles = [0, 45, 90, 135, 180, -135, -90, -45];

% rois 
rois = {'V1', 'V2', 'V3', 'hV4', 'hMTcomplex', 'pMT', 'pMST'};

% data (voxel)
% all

% subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-0397', ...
%     'sub-0442', 'sub-wlsubj121', ...
%     'sub-wlsubj123', 'sub-wlsubj124', 'sub-wlsubj127', 'sub-0395', 'sub-0426', ...
%     'sub-0427', 'sub-0250'};

%% DEFINE the conditions %%
bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
load(fullfile(bidsDir, 'derivatives', sprintf('%sGLM', project), 'hRF_glmsingle', 'allvoxelsBOLDpa.mat'))
load(fullfile(bidsDir, 'derivatives', sprintf('%sGLM', project), 'hRF_glmsingle', 'allparamsBOLDpa.mat'))
% load('/Users/rje257/Desktop/allvoxelsBOLDpa.mat')
% load('/Users/rje257/Desktop/allparamsBOLDpa.mat')

% filter based on a specific eccentricity
eccMin = 4; %1; %0.2;
eccMax = 8; 
rsquaredMin = 0.1; %25;

eccentricity_data = squeeze(allparamsBOLDpa(:,:,:,:,2,:)); % 2 is for eccentricity
conditionExclude = eccentricity_data < eccMin | eccentricity_data >= eccMax;
allvoxelsBOLDpa(conditionExclude) = nan;
rsquared_data = squeeze(allparamsBOLDpa(:,:,:,:,3,:)); % 3 is for r^2
conditionExclude = rsquared_data < rsquaredMin;
allvoxelsBOLDpa(conditionExclude) = nan;

% filter based on subject
subjectStart = 1;
subjectEnd = size(allvoxelsBOLDpa,5);

% filter based on condition
analysisNames = {'Orientation - Baseline', 'Motion - Baseline', 'Motion - Orientation'};
analysisIter = 3;
analysisName = analysisNames{analysisIter};

if analysisIter == 1
    cardArray = 26:27; % all of these variables are motion/orientation
    oblArray = 28:29;
    allArray = [26, 28, 27, 29];
    updownArray = 27;
    rightleftArray = 26;
    upperrightlowerleft = 28;
    upperleftlowerright = 29;
elseif analysisIter == 2
    cardArray = 18:21;
    oblArray = 22:25;
    allArray = [18, 22, 20, 24, 19, 23, 21, 25];
    updownArray = [20, 21];
    rightleftArray = [18,19];
    upperrightlowerleft = [22,23];
    upperleftlowerright = [24,25];
elseif analysisIter == 3
    cardArray = 8:11;
    oblArray = 12:15;
    allArray = [8, 12, 9, 13, 10, 14, 11, 15]; % indices for 0 - 315 in 45 increments
    updownArray = [9,11];
    rightleftArray = [8,10];
    upperrightlowerleft = [12, 14];
    upperleftlowerright = [13, 15];
end

%% plot the BOLD magnitude in the visual field (not accounting for pRF size)
figure
rotationOn = 0; %1; %1; %1;
threshSet = 4; %0.5;
nPoints = 50; %100;

subjectStart = 1;
subjectEnd = size(allvoxelsBOLDpa,5);

if length(subjectStart:subjectEnd)==1
    subName = subjects{subjectStart};
else 
    subName = 'allSubjs';
end

if strcmp(project, 'dg')
    if strcmp(analysisName, 'Orientation - Baseline')
        dirlabels = {'horizontal', 'upright', 'vertical', 'upleft'};
    else
        dirlabels = {'rightwards'; 'upper rightwards'; 'upwards'; 'upper leftwards'; 'leftwards'; ...
            'lower leftwards'; 'downwards'; 'lower rightwards'};
    end
elseif strcmp(project, 'da')
    if strcmp(analysisName, 'Orientation - Baseline')
        dirlabels = {'annulus', 'clockspiral', 'pinwheel', 'cclockspiral'};
    else
        dirlabels = {'clock'; 'clockout'; 'out'; 'cclockout'; 'cclock'; ...
            'cclockin'; 'in'; 'clockin'};
    end
end

for ri=1:length(rois)
    roiNum = ri;
    
    intensity = []; thetaRad = []; rho = []; 
    
    counter = 0; % just to double check # of directions for title
    
    for di=[1,3,5,7] %length(allArray) % 8 directions
        subBETAS = squeeze(allvoxelsBOLDpa(allArray(di), :, roiNum, :, subjectStart:subjectEnd));

        eccen_data = squeeze(allparamsBOLDpa(allArray(di),:,roiNum,:,2,subjectStart:subjectEnd));
        polAng_data = squeeze(allparamsBOLDpa(allArray(di),:,roiNum,:,1,subjectStart:subjectEnd)); % pa ( 0 = UVM; -90 left HM (RH); 90 right HM (LH))
        %polAngMag_data = squeeze(allparamsBOLDpa(:,:,:,:,4,:)); % pa magnitude (NOT CORRECT - DO NOT USE)
        
        % to querry # of vertices for polar angl 1 for this roi, for all
        % subjects (polar angles 1=vertical, 2, 3, 4, 5, 6, 7, 8 goes counterclockwise)
        % squeeze(sum(~isnan(polAng_data(1,:,1:10))))
        
        % allArray(di) is di= 0, ... 315; with 0 as horizontal
        polAng_data_converted = map_theta(polAng_data);
           
        counter = counter+1;
        titlename = dirlabels{di};
        
        if rotationOn
            if di == 1 % rightward motion (rotate counterclockwise 90 deg) ** this one I have to hande the wrapping.
                polAng_data_converted = polAng_data_converted + 90; % values above 360 or below 0 are corrected outside this conditional
            elseif di == 2 % Upper right motion (rotate counterclockwise 45 deg)
                polAng_data_converted = polAng_data_converted + 45;
            elseif di ==3 % Upward motion (no rotation 0 deg)
                polAng_data_converted = polAng_data_converted + 0;
            elseif di == 4 % Upper left motion (rotate clockwise 45 deg)
                polAng_data_converted = polAng_data_converted - 45;
            elseif di == 5 % leftward motion (rotate clockwise 90 deg)
                polAng_data_converted = polAng_data_converted - 90;
            elseif di == 6 % lower leftward (rotate clockwise 135 deg)
                polAng_data_converted = polAng_data_converted - 135;
            elseif di == 7 % downward (rotate clockwise 180 deg)
                polAng_data_converted = polAng_data_converted - 180;
            elseif di == 8 % lower rightward (rotate clockwise 225 deg)
                polAng_data_converted = polAng_data_converted - 225;
            end

            % to fix the "wrapping"
            fixIdx = find(polAng_data_converted>360);
            polAng_data_converted(fixIdx) = polAng_data_converted(fixIdx) -360;
            fixIdx = find(polAng_data_converted<0);
            polAng_data_converted(fixIdx) = 360 + polAng_data_converted(fixIdx);

            if any(polAng_data_converted<0) | any(polAng_data_converted> 360)
                warning('invalid angles')
            end
        end

%         disp('x')

        intensity = [intensity, squeeze(subBETAS(:))];
        thetaRad = [thetaRad, squeeze(deg2rad(polAng_data_converted(:)))];
        rho = [rho, squeeze(eccen_data(:))];
        
    end
    
    % override titlename if more than one direction included
    if counter>1
        titlename = 'all directions';
    end
    
    % after the di loop

    validIdx = find(~isnan(intensity));
    intensity = intensity(validIdx);
    thetaRad = thetaRad(validIdx);
    rho = rho(validIdx);
 
    
    intensity(intensity>threshSet) = threshSet;
    intensity(intensity<(-1*threshSet)) = -1*threshSet;
    
    %intensity = intensity+min(intensity); % normalizing just temporarilty (to compare within plot)
    
    %intensity(intensity<0) = 0;

    subplot(2,4,ri)   
%     %polarplot(thetaRad, rho, 'o');
%     %hold on
% 
    % Loop through the points to plot them with intensity colors
    x = rho .* cos(thetaRad);
    y = rho .* sin(thetaRad);
% 
%     scatter_handle = scatter(x, y, 50, intensity, 'filled'); 
%     scatter_handle.MarkerFaceAlpha = 0.01;
% 
%     % Apply the colormap and color bar
%     colormap(jet);
%     c = colorbar;
%     c.Label.String = 'Intensity';
%     caxis([min(intensity) max(intensity)]);
% 
%     axis equal
%     axis square
% 
%     % Adjust color limits to emphasize positive and negative values
%     caxis([-max(abs(intensity)), max(abs(intensity))]);
% 
%     % Optionally, create a custom colormap to emphasize the contrast
%     % Custom colormap creation
%     n = 128; % Number of color steps
%     neg_colors = [linspace(0, 0, n)', linspace(0, 0.5, n)', linspace(0.5, 1, n)']; % Dark blue to light blue
%     pos_colors = [linspace(0, 1, n)', linspace(0, 0.5, n)', linspace(0, 0, n)']; % Dark red to pink
%     cmap = [neg_colors; flipud(pos_colors)]; % Combine and flip positive colors for smooth transition
% 
%     colormap(cmap);
%     ylim([-15 15])
%     xlim([-15 15])
%     title(rois{roiNum}, 'FontSize',20)


    %z = sin(2*pi*x) .* cos(2*pi*y); % Example data function

    % Create a grid for interpolation
    [Xq, Yq] = meshgrid(linspace(-12, 12, nPoints), linspace(-12, 12, nPoints));

    % Interpolate data onto the grid (averages values at the same location)
    Zq = griddata(x, y, intensity, Xq, Yq, 'cubic'); %'cubic');

    % Convert Cartesian grid to polar coordinates
    %[Theta, Rho] = cart2pol(Xq, Yq);

%     % Plot the data in polar coordinates
%     polaraxes; % Create a polar axes

    % Use polarplot if you want to plot points
    % Alternatively, you can use imagesc or surf for surface plots
    pcolor(Xq, Yq, Zq); % Plot the interpolated data
    shading interp; % Smooth the shading
    colorbar; % Add a color bar
    title(rois{roiNum}, 'FontSize',20);
    
    hold on
    % fill in center with white circle:
    theta = linspace(0, 2*pi, 100);
    r = 4;
    x = r * cos(theta);
    y = r * sin(theta);
    fill(x, y, 'w', 'EdgeColor', 'w') 
    axis equal

    angles_deg = (0:45:315)+22.5;
    angles_rad = deg2rad(angles_deg);
    
    r_start = 4;
    r_end = 8;
    
    % Draw each line
    for i = 1:length(angles_rad)
        x_coords = [r_start, r_end] .* cos(angles_rad(i));
        y_coords = [r_start, r_end] .* sin(angles_rad(i));
        line(x_coords, y_coords, 'Color', 'w', 'LineWidth', 3);
    end
    
    % Optionally plot the annulus again for context
    theta = linspace(0, 2*pi, 200);
    fill([10*cos(theta), 8*cos(fliplr(theta))], [10*sin(theta), 8*sin(fliplr(theta))], 'w', 'EdgeColor', 'w');

    hold on
    plot(a, b, 'k-', 'LineWidth', 1); % 'b-' specifies a blue line
    hold on
    plot(a2, b2, 'k-', 'LineWidth', 1); % 'b-' specifies a blue line
    hold on




    % Create a custom colormap
    n = 256; % Number of colors in the colormap
    blue = [linspace(1, 0, n/2)', linspace(1, 0, n/2)', ones(n/2, 1)]; % Blue to white
    red = [ones(n/2, 1), linspace(0, 1, n/2)', linspace(0, 1, n/2)']; % White to red
    customColormap = [flipud(blue); rgbArray]; %flipud(red)];

    % Apply the custom colormap
    colormap(customColormap);
    %colormap('jet');

    % Ensure the color limits are symmetric around zero
    %caxis([-max(abs(Zq(:))) max(abs(Zq(:)))]);
    caxis([-threshSet threshSet]);
    %caxis([0 threshSet]);
    c = colorbar;
    c.Label.String = 'PSC';

    xlim([-12, 12])
    ylim([-12, 12])

end

f1 = gcf;
f1.Position = [61 345 2252 992];
if analysisIter == 1
    sgtitle(sprintf('%s pRF responses to %s orientation (rotationOn=%i)', subName, titlename, rotationOn))
elseif analysisIter == 2
    sgtitle(sprintf('%s pRF responses to %s motion - baseline (rotationOn=%i)', subName, titlename, rotationOn))
elseif analysisIter == 3
    sgtitle(sprintf('%s pRF responses to %s motion - static (rotationOn=%i)', subName, titlename, rotationOn))
end

%%
figure
for ri=1:length(rois)
    roiNum = ri;
    %subjNum = 1;

    subBETAS_cond1 = squeeze(allvoxelsBOLDpa(cardArray, :, roiNum, :, subjectStart:subjectEnd));
    subBETAS_cond2 = squeeze(allvoxelsBOLDpa(oblArray, :, roiNum, :, subjectStart:subjectEnd));

    subplot(2,4,ri)
    histogram((subBETAS_cond1(:)), 'BinWidth', .01,'FaceAlpha', 0.5, 'FaceColor',[127/255, 191/255, 123/255])
    hold on
    histogram((subBETAS_cond2(:)), 'BinWidth', .01,'FaceAlpha', 0.3, 'FaceColor',[175/255, 141/255, 195/255])
    xlim([-1, 1])
    hold on
    xline(0, 'LineWidth',3)
    title(rois{roiNum}, 'FontSize',20)
end
sgtitle(sprintf('%s Ecc %.10g - %.10g', analysisName, eccMin, eccMax), 'FontSize',20)
f1 = gcf;
f1.Position = [342 589 2184 748];

% radial/tangential

%[90, 45, 0, 315, 270, 225, 180, 135];

figure
for ri=1:length(rois)
    
    roiNum = ri;

    subBETAS_radial1 = squeeze(allvoxelsBOLDpa(updownArray, [1,5], roiNum, :, :)); %up down motion
    subBETAS_tang1 = squeeze(allvoxelsBOLDpa(updownArray, [3,7], roiNum, :, :)); 

    subBETAS_radial2 = squeeze(allvoxelsBOLDpa(rightleftArray, [3,7], roiNum, :, :)); % right left motion
    subBETAS_tang2 = squeeze(allvoxelsBOLDpa(rightleftArray, [1,5], roiNum, :, :)); % right left motion

    subBETAS_radial3 = squeeze(allvoxelsBOLDpa(upperrightlowerleft, [2,6], roiNum, :, :)); % upperright lowerleft motion
    subBETAS_tang3 = squeeze(allvoxelsBOLDpa(upperrightlowerleft, [4,8], roiNum, :, :)); % upperright lowerleft motion

    subBETAS_radial4 = squeeze(allvoxelsBOLDpa(upperleftlowerright, [4,8], roiNum, :, :)); % upperleft lowerright motion
    subBETAS_tang4 = squeeze(allvoxelsBOLDpa(upperleftlowerright, [2,6], roiNum, :, :)); % upperleft lowerright motion

    subBETAS_radial = [subBETAS_radial1(:); subBETAS_radial2(:); subBETAS_radial3(:);subBETAS_radial4(:)];
    subBETAS_tang = [subBETAS_tang1(:); subBETAS_tang2(:); subBETAS_tang3(:);subBETAS_tang4(:)];
    
    
    subplot(2,4,ri)
    histogram((subBETAS_radial(:)), 'BinWidth', .01,'FaceAlpha', 0.5, 'FaceColor',[146 197 222]/255)
    hold on
    histogram((subBETAS_tang(:)), 'BinWidth', .01,'FaceAlpha', 0.2, 'FaceColor',[202 0 32]/255)
    xlim([-1, 1])
    hold on
    xline(0, 'LineWidth',3)
    title(rois{roiNum}, 'FontSize',20)
    %title(sprintf('%s Ecc %.10g - %.10g', rois{roiNum}, eccMin, eccMax), 'FontSize',20)
    
end
sgtitle(sprintf('%s Ecc %.10g - %.10g', analysisName, eccMin, eccMax), 'FontSize',20)
f1 = gcf;
f1.Position = [342 589 2184 748];

%%

num_bins = 50;
concatenated_array = cat(4, subBETAS_cond1, subBETAS_cond2);
[counts, edges] = histcounts(concatenated_array, num_bins);

% Compute the bin centers
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

% Compute the average magnitude within each bin
average_magnitude = zeros(size(counts));
for i = 1:num_bins
    indices = subBETAS_cond1 >= edges(i) & subBETAS_cond1 < edges(i+1);
    average_magnitude(i) = mean(subBETAS_cond1(indices));
end

%[counts2, edges2] = histcounts(subBETAS_cond2, num_bins);

% Compute the bin centers
%bin_centers2 = (edges2(1:end-1) + edges2(2:end)) / 2;

% Compute the average magnitude within each bin
average_magnitude2 = zeros(size(counts));
for i = 1:num_bins
    indices = subBETAS_cond2 >= edges(i) & subBETAS_cond2 < edges(i+1);
    average_magnitude2(i) = mean(subBETAS_cond2(indices));
end


figure
% Plot the average magnitude within each bin
bar(bin_centers, average_magnitude, 'FaceAlpha', 0.3);
hold on
bar(bin_centers, average_magnitude2, 'FaceAlpha', 0.3);
xlabel('Bin Centers');
ylabel('Average Magnitude');
title('Average Magnitude within Each Bin');
hold on
xline(0, 'LineWidth',2)

figure
scatter(bin_centers, average_magnitude-average_magnitude2)
yline(0)

% NOTES:
% orientation
% 90 then 135 then 45 then 0
