%% COMPUTE DERIVATIVES BASED ON POLAR ANGLE

% for DG - this would be both polar cardinal / polar oblique as well as
   % radial / tangential
% for DA - this would be only cartesian cardinal / oblique

function newMatrix = compute_derivativeDirections(medianBOLDpa, projectSettings, varargin)

    % Check project name
    if strcmp(projectSettings.projectName, 'da') || strcmp(projectSettings.projectName, 'dots')
        % Third argument is not required for "da"
        radialvstang = 0;
        asymmetryName = 'derivedCardinalVsDerivedOblique'; % cartesian cardinal

    elseif strcmp(projectSettings.projectName, 'dg')

        % Ensure a third argument is provided for "da"
        if nargin < 3 || isempty(varargin{1})
            error('A third input is required when projectName is "dg".');
        else
            radialvstang = varargin{1};

            if radialvstang == 0
                asymmetryName = 'derivedCardinalVsDerivedOblique'; % polar cardinal
            elseif radialvstang == 1
                asymmetryName = 'radialVsTangential';
            end

        end
    else
        error('Unknown projectName. Expected "da" or "dg".');
    end

    projectName = projectSettings.projectName;
    comparisonName = projectSettings.comparisonName;
    contrasts_dict = projectSettings.contrasts_dict;

    %[proConditions, conConditions, allConditions] = retrieveProConIdx(projectName, comparisonName, radialvstang);

    % I think I always need to assume main cardinal, since radial / tang
    % for dg is defined ad hoc, outside of this function
    [proConditions, conConditions, allConditions] = retrieveProConIdx(projectName, comparisonName, 0);

    % Extract relevant motion direction indices
    motionDirectionIndices = allConditions;
    
    % Polar angles in the index order
    polarAngles = [90, 45, 0, 315, 270, 225, 180, 135];
    %polarAngles = [90, 135, 180, 225, 270, 315, 0, 45];
    
    % Initialize the new matrix:
    % For dg, third dimension is radial, tangential, and other
    % For da, third dimension is vertical, horizontal, and oblique
    newMatrix = zeros(3, size(medianBOLDpa, 2), size(medianBOLDpa, 3), size(medianBOLDpa, 4));
    
    for si = 1:size(medianBOLDpa, 4) % subjects
        for ri = 1:size(medianBOLDpa, 3) % region
    
            for polarIndex = 1:length(polarAngles)
                % Extract data for the current polar angle
                currentPolarData = medianBOLDpa(motionDirectionIndices, polarIndex, ri, si);
                
                provals = [];
                convals = [];
                neutvals = [];
            
                for mi=motionDirectionIndices 
          
                    %%%%%% get the motion direction (in absolute terms)
                        % for dg - this is just the value
                        % for da, later will compute localMotionDirs per
                        % polar angle

                    % retrieve string with motion description
                    contrastValue = contrasts_dict.contrasts(mi).dg_contrast_name;
                    
                    % Extract the part before '_v_'
                    before_v = regexp(contrastValue, '^(.*)_v_', 'tokens', 'once');
                    
                    if ~isempty(before_v) % Extract the numeric part and convert to an integer
                        motionValue = str2double(regexp(before_v{1}, '\d+', 'match', 'once'));
                    else
                        error('No match found in the string');
                    end
                    
                    % Display the result
                    %disp(motionValue); %%%%%%%%

    
                    currentVal = medianBOLDpa(mi, polarIndex, ri, si);
            
                    % check if this is correct - might need to move this out of
                    % the for loop
                    if strcmp(projectName, 'da') || strcmp(projectName, 'dots')
                        UVM_dir = motionValue;
                        localMotionDirs = deriveLocalMotionfromUVM(UVM_dir, polarAngles);
                        currentlocalMotionDir = localMotionDirs(polarIndex);
    
                        if ismember(currentlocalMotionDir, [90, 270])
                            provals = [provals, currentVal]; % this is actually up/down
                        elseif ismember(currentlocalMotionDir, [0, 180])
                            convals = [convals, currentVal]; % this is actually right/left
                        else
                            neutvals = [neutvals, currentVal];
                        end
    
                    elseif strcmp(projectName, 'dg')
                        % Compute the distances between motion directions and polar angles
                        distances = abs(polarAngles(polarIndex) - motionValue);
                        
                        % Identify indices for radial, tangential, and other
                        if distances == 0 || distances == 180 
                            provals = [provals, currentVal]; % radial (in/out)
                        elseif distances == 90 || distances == 270
                            convals = [convals, currentVal]; % tangential (clock/counterclock)
                        else % other
                            neutvals = [neutvals, currentVal];
                        end
                    end
                end
            
                % Compute the averages and store in the new matrix
                newMatrix(1, polarIndex, ri, si) = mean(provals);
                newMatrix(2, polarIndex, ri, si) = mean(convals);
                newMatrix(3, polarIndex, ri, si) = mean(neutvals);
            end
        end
    end

end