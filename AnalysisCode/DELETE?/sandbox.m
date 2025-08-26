%% JUST TO HAVE AS BACKUP -- I THINK I CAN DELETE THIS


    %% COMPUTE DERIVATIVES BASED ON POLAR ANGLE

% Extract relevant motion direction indices
motionDirectionIndices = allConditions; %18:25; %8:15; %26:29; %18:25; %8:15; %8:15; %18:25; %26:29; %

motiondirs = [nan,nan,nan,nan,nan,nan,nan, ...
    0,90,180,270,45,135,225,315, ...
    nan, nan, ...
    0, 180, 90, 270, 45, 225, 135, 315, ...
    0, 90, 45, 135];

% Polar angles in the index order
polarAngles = [90, 45, 0, 315, 270, 225, 180, 135];
%polarAngles = [90, 135, 180, 225, 270, 315, 0, 45];

% Initialize the new matrix for radial, tangential, and other
newMatrix = zeros(3, size(medianBOLDpa, 2), size(medianBOLDpa, 3), size(medianBOLDpa, 4));

for si = 1:size(medianBOLDpa, 4) % subjects
    for ri = 1:size(medianBOLDpa, 3) % region

        for polarIndex = 1:length(polarAngles)
            % Extract data for the current polar angle
            currentPolarData = medianBOLDpa(motionDirectionIndices, polarIndex, ri, si);
            
            radialvals = [];
            tangvals = [];
            othervals = [];
        
            for mi=motionDirectionIndices %8:15 
        
                motionValue = motiondirs(mi);

                currentVal = medianBOLDpa(mi, polarIndex, ri, si);
        
                % check if this is correct - might need to move this out of
                % the for loop
                if strcmp(projectName, 'da')
                    UVM_dir = motionValue;
                    localMotionDirs = deriveLocalMotionfromUVM(UVM_dir, polarAngles);
                    currentlocalMotionDir = localMotionDirs(polarIndex);

                    if ismember(currentlocalMotionDir, [90, 270])
                        radialvals = [radialvals, currentVal]; % this is actually up/down
                    elseif ismember(currentlocalMotionDir, [0, 180])
                        tangvals = [tangvals, currentVal]; % this is actually right/left
                    else
                        othervals = [othervals, currentVal];
                    end

                elseif strcmp(projectName, 'dg')
                    % Compute the distances between motion directions and polar angles
                    distances = abs(polarAngles(polarIndex) - motionValue);
                    
                    % Identify indices for radial, tangential, and other
                    if distances == 0 || distances == 180 % radial
                        radialvals = [radialvals, currentVal];
                    elseif distances == 90 || distances == 270 % tangential
                        tangvals = [tangvals, currentVal];
                    else % other
                        othervals = [othervals, currentVal];
                    end
                end
            end
        
            % Compute the averages and store in the new matrix
            newMatrix(1, polarIndex, ri, si) = mean(radialvals);
            newMatrix(2, polarIndex, ri, si) = mean(tangvals);
            newMatrix(3, polarIndex, ri, si) = mean(othervals);
        end
    end
end