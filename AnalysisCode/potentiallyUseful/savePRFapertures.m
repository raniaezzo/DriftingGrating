% save ret images

experiment = 'benson'; %'benson' for bar or 'dg' for wedgering

if strcmp(experiment, 'benson')
    stimtype = 'RETBAR';
    fileinfo = matfile(sprintf('/Users/rje257/Downloads/apertures/%ssmall.mat', stimtype));
    stim = fileinfo.stim;
    stim = stim./255;
    stim = round(stim);
elseif strcmp(experiment, 'dg')
    stimtype = 'WEGDERING';
    fileinfo = matfile('/Users/rje257/Desktop/transfer_gate/stimfiles.mat');
    stim = fileinfo.stimfiles;
    stim = stim{2};
end


outputFolder = sprintf('/Users/rje257/Downloads/apertures/%ssmallframes', stimtype); % Folder to save images
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder); % Create folder if it doesn't exist
end

for t = 1:size(stim, 3) 
    imshow(stim(:, :, t), []); 
    filename = fullfile(outputFolder, sprintf('frame_%03d.png', t)); % Format filename
    imwrite(stim(:, :, t), filename); % Save as image
end

%%

dutyCycle = compute_duty_cycle_per_cycle(stim);

dutyCycle(dutyCycle==0) = nan;

nanmean(dutyCycle, 'all')

%%

function dutyCycle = compute_duty_cycle_per_cycle(stim)
% Computes the duty cycle for each pixel in an nxnxt stimulus matrix,
% based on ON time within individual cycles.
%
% INPUT:
%   stim - nxnxt binary matrix (1 = ON, 0 = OFF)
%
% OUTPUT:
%   dutyCycle - nxn matrix with the average duty cycle (%) per pixel

    [nx, ny, t] = size(stim); % Get dimensions
    dutyCycle = zeros(nx, ny); % Initialize duty cycle matrix
    
    for x = 1:nx
        for y = 1:ny
            pixelTimeSeries = squeeze(stim(x, y, :)); % Extract time series for pixel (x,y)
            
            % Find all ON indices
            onIndices = find(pixelTimeSeries > 0); 
            
            % Skip if no ON events for this pixel
            if isempty(onIndices)
                dutyCycle(x, y) = 0;
                continue;
            end
            
            % Find cycle start points (first ON event after OFF)
            cycleStartIdx = onIndices([1; find(diff(onIndices) > 1) + 1]); 
            
            % Find cycle durations (time between starts)
            cycleDurations = diff([cycleStartIdx; t]); % Include last cycle duration
            
            % Compute ON time per cycle
            onTimesPerCycle = arrayfun(@(idx, dur) sum(pixelTimeSeries(idx:min(idx + dur - 1, t))), cycleStartIdx, cycleDurations);
            
            % Compute duty cycle per cycle and average over all cycles
            dutyCycle(x, y) = mean(onTimesPerCycle ./ cycleDurations) * 100;
        end
    end
end
