% plot label as image
clear all; clc; close all;

bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';

hemis = {'lh', 'rh'}; %, 'rh'};
roinames = {'pMT', 'pMST', 'hMTcomplex'}; %{'hMTcomplex'};
subjectlist =  {'sub-0427'};

% Define the values and colors
vals = [180, 120, 60, 0]; % Corresponding "val" points
colors = [
    128, 0, 0;    % r, g, b for val=180
    255, 0, 255;  % r, g, b for val=120
    24, 0, 133;   % r, g, b for val=60
    0, 255, 255   % r, g, b for val=0
] / 255; % Normalize to [0, 1]

% Define the number of colors for the colormap
numColors = 255;

% Linearly interpolate the colormap
interpVals = linspace(min(vals), max(vals), numColors);
base = [1 1 1];
rInterp = interp1(vals, colors(:, 1), interpVals, 'linear');
gInterp = interp1(vals, colors(:, 2), interpVals, 'linear');
bInterp = interp1(vals, colors(:, 3), interpVals, 'linear');

% jet map for eccentricity
% eccColormap = jet;
% eccColormap = flipud(eccColormap);
eccvals = [0.2 3.2 6.2 9.2 12.2]; % Corresponding "val" points
%eccvals = [0.2,0.5589,1.5620,4.3654,12.2]; 
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



%%

for ss=1:numel(subjectlist)

    subject=subjectlist{ss};
    fsdir = fullfile(bidsDir, 'derivatives', 'freesurfer', subject);
    prfvistamovdir = retrieve_sesname(fullfile(bidsDir, 'derivatives', 'prfvista_mov_wedgeringbar', subject));


    for ri=1:numel(roinames)
    
        roiname = roinames{ri};
    
        
        for hi=1:numel(hemis)
            roiIndx = [];
    
            hh=hemis{hi};

            % Combine into a colormap
            customColormap = [rInterp', gInterp', bInterp'];
            if strcmp(hh, 'lh')
                customColormap = [base; customColormap];
            elseif strcmp(hh, 'rh')
                customColormap = [flipud(customColormap); base];
            end

            % Load surface points and faces
            [vertices, faces] = freesurfer_read_surf(fullfile(fsdir, 'surf', sprintf('%s.inflated', hh)));

            roi_label = fullfile(fsdir, 'label', 'retinotopy_RE', sprintf('%s.%s_REmanual.label', hh, roiname)); % _REmanual
            roiData = readmatrix(roi_label, 'FileType','text');
            
            roiIndx = [roiIndx; roiData(:,1)+1];


            %faces = faces(roiIndx,:);
            %vertices = vertices(roiIndx,:);

            prfvistamov_eccname = fullfile(prfvistamovdir, sprintf('%s.eccen.mgz', hh));
            eccVol = MRIread(prfvistamov_eccname); eccVol = eccVol.vol;
            %eccVol = eccVol(roiIndx);

            prfvistamov_angname = fullfile(prfvistamovdir, sprintf('%s.angle_adj.mgz', hh));
            angVol = MRIread(prfvistamov_angname); angVol = angVol.vol;

            prfvistamov_vexpname = fullfile(prfvistamovdir, sprintf('%s.vexpl.mgz', hh));
            veVol = MRIread(prfvistamov_vexpname); veVol = veVol.vol;

            % Extract RAS coordinates and label values
            x = roiData(:, 2); % X-coordinates in RAS
            y = roiData(:, 3); % Y-coordinates in RAS
            z = roiData(:, 4); % Z-coordinates in RAS (optional for 2D)
            labelValues = roiData(:, 5); % Label values

%             if ri==3
%                 figure
%                 % Visualize the surface with eccentricity values
%                 trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), eccVol', 'EdgeColor', 'none');
%                 colormap(jet);
%                 colorbar;
%                 caxis([0, 12]); 
%                 title('Eccentricity on Flattened Surface');
%                 axis equal;
%             end
                
            % Create a mask for the label
            mask = ismember(1:size(vertices, 1), roiIndx);
            
            % Filter vertices
            verticesROI = vertices(mask, :);
            
            % Find faces within the ROI
            facesROI = faces(all(mask(faces), 2), :); % Faces where all vertices are in the label
            
            % Filter the eccentricity values
            eccVolROI = eccVol(mask);
            angVolROI = angVol(mask);
            veVolROI = veVol(mask);
            
            % Use PCA to define the principal plane
            [coeff, ~, ~] = pca(verticesROI);
            flattenedCoords = verticesROI * coeff(:, 1:2); % Project onto the first two principal components
            
            % Define the 2D grid resolution
            gridSize = 100; % Adjust as needed
            xRange = linspace(min(flattenedCoords(:, 1)), max(flattenedCoords(:, 1)), gridSize);
            yRange = linspace(min(flattenedCoords(:, 2)), max(flattenedCoords(:, 2)), gridSize);
            [xGrid, yGrid] = meshgrid(xRange, yRange);
            
            % Interpolate the eccentricity values onto the 2D grid
            eccGrid{hi,ri} = griddata(flattenedCoords(:, 1), flattenedCoords(:, 2), eccVolROI, xGrid, yGrid, 'cubic');
            eccGrid{hi,ri}(eccGrid{hi,ri} == 0) = NaN;
            
            angGrid{hi,ri} = griddata(flattenedCoords(:, 1), flattenedCoords(:, 2), angVolROI, xGrid, yGrid, 'cubic');
            angGrid{hi,ri}(angGrid{hi,ri} == 0) = NaN;

            veGrid{hi,ri} = griddata(flattenedCoords(:, 1), flattenedCoords(:, 2), veVolROI, xGrid, yGrid, 'cubic');
            veGrid{hi,ri}(veGrid{hi,ri} == 0) = NaN;
            
            if ri==3
                figure;
                subplot(1,2,1)
                h = imagesc(xRange, yRange, eccGrid{hi, ri}); % Display the flattened matrix
                ax2 = gca;
                colormap(ax2, eccColormap);
                colorbar;
                %caxis([0, 30]);
                caxis([0, 8]); 
                %caxis([1, 8]); %caxis([0, 30]);
                xlabel('Flattened X');
                ylabel('Flattened Y');
                title('Flattened Region with Eccentricity Values');
                axis equal;
                xticks([])
                yticks([])
                set(gca, 'Color', 'w'); % Set background color to white
                alphaData = double(~isnan(eccGrid{hi, ri})); % Create AlphaData: 1 for valid, 0 for NaN
                set(h, 'AlphaData', alphaData); 
                hold off
                
                
                subplot(1,2,2)
                h = imagesc(xRange, yRange, angGrid{hi, ri}); % Display the flattened matrix
                ax1 = gca; % Get the handle for the first subplot
                colormap(ax1, customColormap); 
                colorbar
                if strcmp(hh, 'lh')
                    clim(ax1, [0 180]); 
                elseif strcmp(hh, 'rh')
                    clim(ax1, [-180 0]); 
                end
                xlabel('Flattened X');
                ylabel('Flattened Y');
                title('Flattened Region with Angle Values');
                axis equal;
                xticks([])
                yticks([])
                set(gca, 'Color', 'w'); % Set background color to white
                alphaData = double(~isnan(angGrid{hi, ri})); % Create AlphaData: 1 for valid, 0 for NaN
                set(h, 'AlphaData', alphaData); 
                sgtitle(sprintf('%s %s',hh, subject))
                f1 = gcf;
                f1.Position = [1 391 1329 946];

                figure
                h = imagesc(xRange, yRange, veGrid{hi, ri}); % Display the flattened matrix
                ax2 = gca;
                colormap('jet');
                colorbar;
                %caxis([0, 30]);
                caxis([0, 1]); 
                %caxis([1, 8]); %caxis([0, 30]);
                xlabel('Flattened X');
                ylabel('Flattened Y');
                title(sprintf('%s Flattened Region with R^2 Values', hh));
                axis equal;
                xticks([])
                yticks([])
                set(gca, 'Color', 'w'); % Set background color to white
                alphaData = double(~isnan(veGrid{hi, ri})); % Create AlphaData: 1 for valid, 0 for NaN
                set(h, 'AlphaData', alphaData); 
                hold off
            end
         end
    end
end

%% now find edge between MT and MST

figure
imagesc(angGrid{1,1})
figure
imagesc(angGrid{1,2})

%% now smooth

if strcmp(hh, 'lh')
    angGrid = angGrid{1,3};
    eccGrid = eccGrid{1,3};
elseif strcmp(hh, 'rh')
    angGrid = angGrid{2,3};
    eccGrid = eccGrid{2,3};
end

% Define the size of the moving average filter
filterSize = 5; % 5x5 window

% Create the filter
kernel = ones(filterSize) / (filterSize^2);

% Apply the moving average filter
smoothedData = imfilter(angGrid, kernel, 'replicate');

%%
% Define the ROI as non-NaN values
%ROI_mask = ~isnan(smoothedData);
ROI_mask = ~isnan(angGrid);

% Find the boundary of the ROI
boundaryStruct = bwboundaries(ROI_mask); % Extract boundary points
boundaryPoints = boundaryStruct{1};     % Boundary points as [row, col]

% Convert to X, Y coordinates
boundaryX = boundaryPoints(:, 2); % Columns correspond to X
boundaryY = boundaryPoints(:, 1); % Rows correspond to Y

% Calculate the width of the ROI
ROI_width = sqrt((max(boundaryX) - min(boundaryX))^2 + (max(boundaryY) - min(boundaryY))^2);
minLength = ROI_width / 3; % Minimum length constraint

% Define 100 evenly spaced points along the boundary
numPoints = 100;
boundaryStep = round(linspace(1, length(boundaryX), numPoints));
boundaryXSampled = boundaryX(boundaryStep);
boundaryYSampled = boundaryY(boundaryStep);

% Initialize variables for finding the optimal polynomial
bestDeviation = Inf;
bestStartIdx = 0;
bestEndIdx = 0;
bestPolyCoeffs = [];

% Loop over all combinations of start and end points
for startIdx = 1:numPoints
    for endIdx = 1:numPoints
        if startIdx == endIdx
            continue; % Skip if start and end points are the same
        end

        % Get start and end points
        startX = boundaryXSampled(startIdx);
        startY = boundaryYSampled(startIdx);
        endX = boundaryXSampled(endIdx);
        endY = boundaryYSampled(endIdx);

        % Calculate the length of the candidate line
        lineLength = sqrt((endX - startX)^2 + (endY - startY)^2);
        if lineLength < minLength
            continue; % Skip if the line is too short
        end

        % Generate raw line between start and end points
        lineX = linspace(startX, endX, 100);
        lineY = linspace(startY, endY, 100);

        % Interpolate eccentricity values along the line
        lineValues = interp2(1:size(smoothedData, 2), 1:size(smoothedData, 1), smoothedData, lineX, lineY, 'linear', NaN);

        % Remove NaNs for polynomial fitting
        validIndices = ~isnan(lineValues);
        if sum(validIndices) > 2 % Ensure sufficient data for polynomial fit
            lineValues = lineValues(validIndices); % Filter valid values
            fitX = 1:length(lineValues); % Independent variable (e.g., position along the line)

            % Fit a 2nd-order polynomial
            coeffs = polyfit(fitX, lineValues, 2); % Fit a 2nd-order polynomial
            fitY = polyval(coeffs, fitX);          % Evaluate the polynomial

            % Compute the deviation from a zero-mean solution
            meanDeviation = abs(mean(fitY)); % Deviation from zero-mean solution

            % Check if this is the best fit
            if meanDeviation < bestDeviation
                bestDeviation = meanDeviation;
                bestStartIdx = startIdx;
                bestEndIdx = endIdx;
                bestPolyCoeffs = coeffs;
            end
        end
    end
end

% Output the best start and end indices
fprintf('Best Start Index: %d, Best End Index: %d\n', bestStartIdx, bestEndIdx);
fprintf('Best Polynomial Coefficients: [%.4f, %.4f, %.4f]\n', bestPolyCoeffs);

% Plot the best polynomial on top of the ROI
figure;
subplot(1,2,1)
imagesc(smoothedData);
ax1 = gca; % Get the handle for the first subplot
colormap(ax1, customColormap); 
if strcmp(hh, 'lh')
    clim(ax1, [0 180]); 
elseif strcmp(hh, 'rh')
    clim(ax1, [-180 0]); 
end
colorbar;
hold on;
plot(boundaryX, boundaryY, 'w-', 'LineWidth', 1); % Plot the ROI boundary
plot(boundaryXSampled, boundaryYSampled, 'k.', 'MarkerSize', 15); % Plot the sampled boundary points
xlabel('X');
ylabel('Y');
title('Optimal Polynomial Fit with Length Constraint');
axis equal;

% Generate the best polynomial line
bestStartX = boundaryXSampled(bestStartIdx);
bestStartY = boundaryYSampled(bestStartIdx);
bestEndX = boundaryXSampled(bestEndIdx);
bestEndY = boundaryYSampled(bestEndIdx);

% Generate the line from best start and end points
bestLineX = linspace(bestStartX, bestEndX, 100);
bestLineY = linspace(bestStartY, bestEndY, 100);

% Evaluate the polynomial for the best line
fitYBest = polyval(bestPolyCoeffs, 1:length(bestLineX));

% Plot the polynomial fit
plot(bestLineX, bestLineY, 'r-', 'LineWidth', 2); % Plot the fitted polynomial

hold off
subplot(1,2,2)
imagesc(eccGrid);
hold on;
plot(boundaryX, boundaryY, 'w-', 'LineWidth', 1); % Plot the ROI boundary
plot(boundaryXSampled, boundaryYSampled, 'k.', 'MarkerSize', 15); % Plot the sampled boundary points
ax2 = gca; % Get the handle for the first subplot
colormap(ax2, jet); 
clim(ax2, [0, 30]); 
colorbar;
axis equal;
