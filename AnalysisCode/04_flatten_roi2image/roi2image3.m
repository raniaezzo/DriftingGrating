%clear all; close all; clc

% Path to the flattened patch file
hemi = 'lh';
subjectname = 'sub-wlsubj127';

%bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids';
bidsDir = '/Volumes/server/Projects/Project_dg/data_bids/';

flatFile = fullfile(bidsDir, sprintf('/derivatives/freesurfer/%s/surf/%s.hMTcomplex.patch.flat', subjectname, hemi));

% Read the flattened patch file
patch = read_patch(flatFile);

patch.y_new = patch.x; patch.x_new = patch.y;
patch.x = patch.x_new; patch.y = patch.y_new; 

% Scatter plot of the patch
figure;
scatter(patch.x, patch.y, 10, 'filled'); % Adjust marker size as needed
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Flattened Patch Visualization');

%%

% Triangulate the patch vertices (2D as z is likely 0)
tri = delaunay(patch.x, patch.y);

% Plot the reconstructed surface
figure;
trisurf(tri, patch.x, patch.y, patch.z, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Flattened Patch Reconstructed Surface');

%%

prfvistaFile = fullfile(bidsDir, sprintf('/derivatives/prfvista_mov/%s/ses-nyu3t01/%s.angle_adj.mgz', subjectname, hemi));

prfvistamov_angname = fullfile(prfvistaFile);
angVol = MRIread(prfvistamov_angname); angVol = angVol.vol;

prfvistaFile = fullfile(bidsDir, sprintf('/derivatives/prfvista_mov/%s/ses-nyu3t01/%s.eccen.mgz', subjectname, hemi));

prfvistamov_eccname = fullfile(prfvistaFile);
eccVol = MRIread(prfvistamov_eccname); eccVol = eccVol.vol;

% Define the values and colors
vals = [180, 120, 60, 0]; % Corresponding "val" points
colors = [
    128, 0, 0;    % r, g, b for val=180
    255, 0, 255;  % r, g, b for val=120
    24, 0, 133;   % r, g, b for val=60
    0, 255, 255   % r, g, b for val=0
] / 255; % Normalize to [0, 1]

% Define the number of colors for the colormap
numColors = 256;

% Linearly interpolate the colormap
interpVals = linspace(min(vals), max(vals), numColors);
base = []; %[1 1 1];
rInterp = interp1(vals, colors(:, 1), interpVals, 'linear');
gInterp = interp1(vals, colors(:, 2), interpVals, 'linear');
bInterp = interp1(vals, colors(:, 3), interpVals, 'linear');

% Combine into a colormap
customColormap = [rInterp', gInterp', bInterp'];
if strcmp(hemi, 'lh')
    customColormap = [base; customColormap];
elseif strcmp(hemi, 'rh')
    customColormap = [base; flipud(customColormap)];
end


% jet map for eccentricity
%eccColormap = jet;
%eccColormap = flipud(eccColormap);
eccvals = [0.2 3.2 6.2 9.2 12.2]; % Corresponding "val" points
ecccolors = [
    255, 0, 0; % red
    255, 255, 0;  
    0, 255, 0;  
    0, 255,255;  
    0, 0, 255;   % blue
] / 255; % Normalize to [0, 1]

ecc_interpVals = linspace(min(eccvals), max(eccvals), numColors);
eccrInterp = interp1(eccvals, ecccolors(:, 1), ecc_interpVals, 'linear');
eccgInterp = interp1(eccvals, ecccolors(:, 2), ecc_interpVals, 'linear');
eccbInterp = interp1(eccvals, ecccolors(:, 3), ecc_interpVals, 'linear');

eccColormap = [eccrInterp', eccgInterp', eccbInterp'];

% Assume 'functionalData' contains vertex-wise data for the original surface
overlayData_ang = angVol(patch.ind+1); % Extract data for the patch vertices
overlayData_ecc = eccVol(patch.ind+1); % Extract data for the patch vertices

% Visualize overlay data on the patch
figure;
scatter(patch.x, patch.y, 100, overlayData_ang, 'filled');
axis equal;
colormap(customColormap);
ax1 = gca; 
if strcmp(hemi, 'lh')
    clim(ax1, [0 180]); 
elseif strcmp(hemi, 'rh')
    clim(ax1, [-180 0]); 
end
colorbar;
xlabel('X');
ylabel('Y');
title(sprintf('%s %s', subjectname, hemi));

%%
% Define grid resolution
xGrid = linspace(min(patch.x), max(patch.x), 100);
yGrid = linspace(min(patch.y), max(patch.y), 100);

% Create the grid
[XGrid, YGrid] = meshgrid(xGrid, yGrid);

% Interpolate overlay data onto the grid
gridData_ang = griddata(patch.x, patch.y, overlayData_ang, XGrid, YGrid, 'linear');
gridData_ecc = griddata(patch.x, patch.y, overlayData_ecc, XGrid, YGrid, 'linear');

% % Visualize the interpolated data
% figure;
% subplot(1,2,1)
% pcolor(XGrid, YGrid, gridData_ang);
% shading interp; % Smooth interpolation between grid points
% axis equal;
% xlabel('X');
% ylabel('Y');
% ax2 = gca; % Get the handle for the first subplot
% colormap(ax2, customColormap); 
% ax1 = gca; 
% if strcmp(hemi, 'lh')
%     clim(ax1, [0 180]); 
% elseif strcmp(hemi, 'rh')
%     clim(ax1, [-180 0]); 
% end
% colorbar;
% title(sprintf('%s %s polar angle', subjectname, hemi));
% 
% subplot(1,2,2)
% pcolor(XGrid, YGrid, gridData_ecc);
% shading interp; % Smooth interpolation between grid points
% axis equal;
% xlabel('X');
% ylabel('Y');
% ax2 = gca; % Get the handle for the first subplot
% colormap(ax2, jet); 
% ax1 = gca; 
% clim(ax2, [0.2, 8]); 
% colorbar;
% title(sprintf('%s %s eccentricity', subjectname, hemi));
% 
% 
% % %% another way to do it
% subplot(1,2,2)
% imagesc(xGrid, yGrid, gridData_ecc);
% alphaData = ~isnan(gridData_ecc); % AlphaData: 1 for non-NaN, 0 for NaN
% set(gca, 'Color', [1, 1, 1]); % Set axes background to white
% set(get(gca, 'Children'), 'AlphaData', alphaData); % Apply transparency
% axis equal;
% xlabel('X');
% ylabel('Y');
% ax2 = gca; % Get the handle for the first subplot
% colormap(ax2, jet); 
% clim(ax2, [0.2, 8]); 
% colorbar;
% title('2D Patch Visualization Using imagesc');

%%

% Next, load in the label (MT, MST) - 
% 1- project onto the patch 

hh = hemi;

fsdir = fullfile(bidsDir, sprintf('/derivatives/freesurfer/%s/', subjectname));

[MT_filteredX, MT_filteredY, MT_vals] = retrieve_roiIdx(fsdir, hh, 'pMT', patch);
[MST_filteredX, MST_filteredY, MST_vals] = retrieve_roiIdx(fsdir, hh, 'pMST', patch);
[MTcomplex_filteredX, MTcomplex_filteredY, MTcomplex_vals] = retrieve_roiIdx(fsdir, hh, 'hMTcomplex', patch);

MTData_ang = griddata(MT_filteredX, MT_filteredY, angVol(MT_vals), XGrid, YGrid, 'linear');
MTData_ecc = griddata(MT_filteredX, MT_filteredY, eccVol(MT_vals), XGrid, YGrid, 'linear');
MSTData_ang = griddata(MST_filteredX, MST_filteredY, angVol(MST_vals), XGrid, YGrid, 'linear');
MSTData_ecc = griddata(MST_filteredX, MST_filteredY, eccVol(MST_vals), XGrid, YGrid, 'linear');
MTcomplexData_ang = griddata(MTcomplex_filteredX, MTcomplex_filteredY, angVol(MTcomplex_vals), XGrid, YGrid, 'linear');
MTcomplexData_ecc = griddata(MTcomplex_filteredX, MTcomplex_filteredY, eccVol(MTcomplex_vals), XGrid, YGrid, 'linear');

% Visualize the interpolated data
figure;
subplot(1,2,1)
pcolor(XGrid, YGrid, MTcomplexData_ang);
shading interp; % Smooth interpolation between grid points
axis equal;
xlabel('X');
ylabel('Y');
ax2 = gca; % Get the handle for the first subplot
colormap(ax2, customColormap); 
ax1 = gca; 
if strcmp(hemi, 'lh')
    clim(ax1, [0 180]); 
elseif strcmp(hemi, 'rh')
    clim(ax1, [-180 0]); 
end
colorbar;
title(sprintf('%s %s polar angle', subjectname, hemi));

subplot(1,2,2)
pcolor(XGrid, YGrid, MTcomplexData_ecc);
shading interp; % Smooth interpolation between grid points
axis equal;
xlabel('X');
ylabel('Y');
ax2 = gca; % Get the handle for the first subplot
colormap(ax2, eccColormap); 
ax1 = gca; 
clim(ax2, [0.2, 8]); 
colorbar;
title(sprintf('%s %s eccentricity', subjectname, hemi));

sgtitle('MTcomplex')

%%
% 2- isolate the boundary of MT/MST (vertical meridian)

figure
scatter(MST_filteredX, MST_filteredY)
hold on
scatter(MT_filteredX, MT_filteredY, 'r')

[slope, intercept] = findSeparationLine(MT_filteredX', MT_filteredY', MST_filteredX', MST_filteredY');

% xVals = linspace(min([MT_filteredX'; MST_filteredX']), max([MT_filteredX'; MST_filteredX']), 100);
% yVals = slope * xVals + intercept;

% Compute y-values range for plotting
minY = min([MT_filteredY, MST_filteredY]);
maxY = max([MT_filteredY, MST_filteredY]);
yVals = linspace(minY, maxY, 100);

% Compute x-values from y-values using the line equation
xVals = (yVals - intercept) / slope;

plot(xVals, yVals, 'k-', 'LineWidth', 2);


%%
% 3- rotate entire image by that vector

[XGridRot_pa, YGridRot_pa, rotatedData_pa] = rotateMap(MTcomplexData_ang, slope);
[XGridRot_ecc, YGridRot_ecc, rotatedData_ecc] = rotateMap(MTcomplexData_ecc, slope);

% Plot the rotated data
figure;
subplot(1,2,1)
pcolor(XGridRot_pa, YGridRot_pa, rotatedData_pa);
shading interp; % Smooth visualization
ax2 = gca; % Get the handle for the first subplot
colormap(ax2, customColormap); 
ax1 = gca; 
if strcmp(hemi, 'lh')
    clim(ax1, [0 180]); 
elseif strcmp(hemi, 'rh')
    clim(ax1, [-180 0]); 
end
axis equal;
xlabel('X');
ylabel('Y');
title('Rotated Data Plot');

subplot(1,2,2)
pcolor(XGridRot_ecc, YGridRot_ecc, rotatedData_ecc);
shading interp; % Smooth visualization
ax2 = gca; % Get the handle for the first subplot
colormap(ax2, eccColormap); 
ax1 = gca; 
clim(ax2, [0.2, 8]); 
axis equal;
xlabel('X');
ylabel('Y');
title('Rotated Data Plot');

%%

% Example data
[XGridRot, YGridRot] = meshgrid(1:120, 1:120);

% Apply the transformation and cropping function
croppedData_pa = transformAndCropForImagesc(XGridRot, YGridRot, rotatedData_pa);
croppedData_pa = flipud(croppedData_pa);
croppedData_ecc = transformAndCropForImagesc(XGridRot, YGridRot, rotatedData_ecc);
croppedData_ecc = flipud(croppedData_ecc);

% Plot using imagesc
figure;
subplot(1,2,1)
imagesc(croppedData_pa);
ax1 = gca; % Get the handle for the first subplot
colormap(ax1, customColormap); 
alphaData = ~isnan(croppedData_pa); % AlphaData: 1 for non-NaN, 0 for NaN
set(ax1, 'Color', [1, 1, 1]); % Set axes background to white
set(get(ax1, 'Children'), 'AlphaData', alphaData); % Apply transparency
if strcmp(hemi, 'lh')
    clim(ax1, [0 180]); 
elseif strcmp(hemi, 'rh')
    clim(ax1, [-180 0]); 
end
axis equal; axis square;
title('Polar angle rotated');

subplot(1,2,2)
imagesc(croppedData_ecc);
ax2 = gca; % Get the handle for the first subplot
colormap(ax2, eccColormap); 
alphaData = ~isnan(croppedData_ecc); % AlphaData: 1 for non-NaN, 0 for NaN
set(ax2, 'Color', [1, 1, 1]); % Set axes background to white
set(get(ax2, 'Children'), 'AlphaData', alphaData); % Apply transparency
clim(ax2, [0.2, 8]); 
axis equal; axis square;
title('Eccentricity rotated');

%%

data = croppedData_pa;
% Stretch data into a 50x50 square
stretchedData = stretchMatrixToSquare(data, 50);

% Plot the results
figure;

subplot(2, 2, 1);
imagesc(data);
axis equal;
ax1 = gca; % Get the handle for the first subplot
colormap(ax1, customColormap); 
alphaData = ~isnan(data); % AlphaData: 1 for non-NaN, 0 for NaN
set(ax1, 'Color', [1, 1, 1]); % Set axes background to white
set(get(ax1, 'Children'), 'AlphaData', alphaData); % Apply transparency
if strcmp(hemi, 'lh')
    clim(ax1, [0 180]); 
elseif strcmp(hemi, 'rh')
    clim(ax1, [-180 0]); 
end
title('Original Data');

subplot(2, 2, 2);
imagesc(stretchedData);
axis equal;
ax2 = gca; % Get the handle for the first subplot
colormap(ax2, customColormap); 
alphaData = ~isnan(stretchedData); % AlphaData: 1 for non-NaN, 0 for NaN
set(ax2, 'Color', [1, 1, 1]); % Set axes background to white
set(get(ax2, 'Children'), 'AlphaData', alphaData); % Apply transparency
if strcmp(hemi, 'lh')
    clim(ax2, [0 180]); 
elseif strcmp(hemi, 'rh')
    clim(ax2, [-180 0]); 
end
title('Stretched Data (50x50 Square)');
h1 = colorbar('Position', [0.92, 0.55, 0.02, 0.3]);


data = croppedData_ecc;
% Stretch data into a 50x50 square
stretchedData = stretchMatrixToSquare(data, 50);

subplot(2, 2, 3);
imagesc(data);
axis equal;
ax3 = gca; % Get the handle for the first subplot
colormap(ax3, eccColormap); 
alphaData = ~isnan(data); % AlphaData: 1 for non-NaN, 0 for NaN
set(ax3, 'Color', [1, 1, 1]); % Set axes background to white
set(get(ax3, 'Children'), 'AlphaData', alphaData); % Apply transparency
clim(ax3, [1, 8]); 
title('Original Data');

subplot(2, 2, 4);
imagesc(stretchedData);
axis equal;
ax4 = gca; % Get the handle for the first subplot
colormap(ax4, eccColormap); 
alphaData = ~isnan(stretchedData); % AlphaData: 1 for non-NaN, 0 for NaN
set(ax4, 'Color', [1, 1, 1]); % Set axes background to white
set(get(ax4, 'Children'), 'AlphaData', alphaData); % Apply transparency
clim(ax4, [0.2, 8]); 
title('Stretched Data (50x50 Square)');
h2 = colorbar('Position', [0.92, 0.1, 0.02, 0.3]);

%%

function [filteredX, filteredY, vals] = retrieve_roiIdx(fsdir, hh, roiname, patch)
    roi_label = fullfile(fsdir, 'label', 'retinotopy_RE', sprintf('%s.%s_REmanual.label', hh, roiname)); % _REmanual
    roiData = readmatrix(roi_label, 'FileType','text');
    
    roiIndx = roiData(:,1)+1;
    
    % create 
    [~, idx] = ismember(patch.ind',roiIndx);
    validIdx = idx > 0;  % this is smaller than the actual ROI - why?
    filteredX = patch.x(validIdx);
    filteredY = patch.y(validIdx);
    
    % Interpolate overlay data onto the grid
    vals = patch.ind(validIdx);
end

function [slope, intercept] = findSeparationLine(X1, Y1, X2, Y2)
    % Combine the data
    X = [X1; X2];
    Y = [Y1; Y2];
    labels = [ones(size(X1)); -ones(size(X2))]; % Group 1 as +1, Group 2 as -1

    % Compute group means
    meanX1 = mean(X1);
    meanY1 = mean(Y1);
    meanX2 = mean(X2);
    meanY2 = mean(Y2);

    % Center the data
    X1Centered = X1 - meanX1;
    Y1Centered = Y1 - meanY1;
    X2Centered = X2 - meanX2;
    Y2Centered = Y2 - meanY2;

    % Compute within-class scatter matrices
    S1 = cov(X1Centered, Y1Centered);
    S2 = cov(X2Centered, Y2Centered);
    Sw = S1 + S2; % Within-class scatter matrix

    % Compute the direction vector for the line
    meanDiff = [meanX1 - meanX2; meanY1 - meanY2];
    w = Sw \ meanDiff; % Linear discriminant direction

    % Derive the slope and intercept of the decision boundary
    slope = -w(1) / w(2);
    intercept = (meanY1 + meanY2) / 2 - slope * (meanX1 + meanX2) / 2;
end

function [XGridRot, YGridRot, rotatedData] = rotateMap(MTcomplexData_ecc, slope)
    [XGrid, YGrid] = meshgrid(1:size(MTcomplexData_ecc, 2), 1:size(MTcomplexData_ecc, 1));
    angleDegrees = atan(slope) * (180 / pi) - 90; % Convert radians to degrees
    
    % Pad the matrix with NaNs to avoid cutting off points
    paddedData = padarray(MTcomplexData_ecc, [10, 10], NaN, 'both');
    rotatedData = imrotate(paddedData, angleDegrees, 'bilinear', 'crop');
    rotatedData(rotatedData == 0) = NaN;
    [XGridRot, YGridRot] = meshgrid(1:size(rotatedData, 2), 1:size(rotatedData, 1));
end

function croppedData = transformAndCropForImagesc(XGridRot, YGridRot, rotatedData)
    % Define a regular grid based on the range of X and Y
    xRange = linspace(min(XGridRot(:)), max(XGridRot(:)), size(rotatedData, 2));
    yRange = linspace(min(YGridRot(:)), max(YGridRot(:)), size(rotatedData, 1));
    [XGridReg, YGridReg] = meshgrid(xRange, yRange);

    % Interpolate the data onto the regular grid
    transformedData = griddata(XGridRot, YGridRot, rotatedData, XGridReg, YGridReg, 'linear');

    % Explicitly replace zeros with NaN where no valid interpolation occurred
    % This ensures zeros aren't misinterpreted as valid data
    transformedData(isnan(rotatedData)) = NaN;

    % Find rows and columns that are fully NaN
    validRows = any(~isnan(transformedData), 2); % Rows with at least one non-NaN
    validCols = any(~isnan(transformedData), 1); % Columns with at least one non-NaN

    % Crop the transformed data to remove NaN-only rows and columns
    croppedData = transformedData(validRows, validCols);
end

function stretchedData = stretchMatrixToSquare(originalData, squareSize)
    % Input:
    % originalData - Input matrix with arbitrary shape (including NaNs)
    % squareSize - Size of the output square (e.g., 50 for 50x50)

    % Get the size of the original data
    [rows, cols] = size(originalData);

    % Create the original coordinate grid
    [X, Y] = meshgrid(1:cols, 1:rows);

    % Flatten the matrix into a set of points
    x = X(:); % All X coordinates
    y = Y(:); % All Y coordinates
    values = originalData(:); % All corresponding values

    % Normalize the coordinates to fit within [0, 1]
    xNorm = (x - min(x)) / (max(x) - min(x));
    yNorm = (y - min(y)) / (max(y) - min(y));

    % Define a uniform square grid
    [Xq, Yq] = meshgrid(linspace(0, 1, squareSize), linspace(0, 1, squareSize));

    % Interpolate onto the square grid
    stretchedData = griddata(xNorm, yNorm, values, Xq, Yq, 'linear');

    % Fill any remaining NaNs (optional)
    stretchedData = fillmissing(stretchedData, 'nearest');
end
