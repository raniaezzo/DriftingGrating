function newMatrix = compute_derivativeEachAbsDirection(medianBOLDpa, projectSettings)

    % Check project name
    if strcmp(projectSettings.projectName, 'da')
        
        projectName = projectSettings.projectName;
        comparisonName = projectSettings.comparisonName;
        contrasts_dict = projectSettings.contrasts_dict;
    
        %[proConditions, conConditions, allConditions] = retrieveProConIdx(projectName, comparisonName, radialvstang);
    
        % I think I always need to assume main cardinal, since radial / tang
        % for dg is defined ad hoc, outside of this function
        allConditions = retrieveAbsDirIdx(projectName, comparisonName);
    
        % Extract relevant motion direction indices
        motionDirectionIndices = allConditions;
        nVals = size(motionDirectionIndices,2);
        
        % Polar angles in the index order
        polarAngles = [90, 45, 0, 315, 270, 225, 180, 135];
        %polarAngles = [90, 135, 180, 225, 270, 315, 0, 45];
        
        % Initialize the new matrix:
        % For dg, third dimension is radial, tangential, and other
        % For da, third dimension is vertical, horizontal, and oblique
        newMatrix = zeros(nVals, size(medianBOLDpa, 2), size(medianBOLDpa, 3), size(medianBOLDpa, 4));
        
        for si = 1:size(medianBOLDpa, 4) % subjects
            for ri = 1:size(medianBOLDpa, 3) % region
        
                for polarIndex = 1:length(polarAngles)
                    % Extract data for the current polar angle
                    currentPolarData = medianBOLDpa(motionDirectionIndices, polarIndex, ri, si);
                    
                    if strcmp(comparisonName, 'orientation_minus_baseline')
                        horizontal = [];
                        upperright = [];
                        vertical = [];
                        upperleft = [];
                    else
                        rightwards = [];
                        upperrightwards = [];
                        upwards = [];
                        upperleftwards = [];
                        leftwards = [];
                        lowerleftwards = [];
                        downwards = [];
                        lowerrightwards = [];
                    end
                
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
                
                        UVM_dir = motionValue;
                        localMotionDirs = deriveLocalMotionfromUVM(UVM_dir, polarAngles);
                        currentlocalMotionDir = localMotionDirs(polarIndex);
    
                        % here orientations include 0-316 due to function
                        % conversion (deriveLocalMotionfromUVM)
                        if strcmp(comparisonName, 'orientation_minus_baseline')
                            if ismember(currentlocalMotionDir, [0,180])
                                horizontal = [horizontal, currentVal];
                            elseif ismember(currentlocalMotionDir, [45, 225])
                                upperright = [upperright, currentVal];
                            elseif ismember(currentlocalMotionDir, [90, 270])
                                vertical = [vertical, currentVal];
                            elseif ismember(currentlocalMotionDir, [135, 315])
                                upperleft = [upperleft, currentVal];
                            else
                                warning('does not meet a condition')
                            end
                        else
                            if ismember(currentlocalMotionDir, 0)
                                rightwards = [rightwards, currentVal];
                            elseif ismember(currentlocalMotionDir, 45)
                                upperrightwards = [upperrightwards, currentVal]; 
                            elseif ismember(currentlocalMotionDir, 90)
                                upwards = [upwards, currentVal];
                            elseif ismember(currentlocalMotionDir, 135)
                                upperleftwards = [upperleftwards, currentVal]; 
                            elseif ismember(currentlocalMotionDir, 180)
                                leftwards = [leftwards, currentVal]; 
                            elseif ismember(currentlocalMotionDir, 225)
                                lowerleftwards = [lowerleftwards, currentVal];
                            elseif ismember(currentlocalMotionDir, 270)
                                downwards = [downwards, currentVal];
                            elseif ismember(currentlocalMotionDir, 315)
                                lowerrightwards = [lowerrightwards, currentVal];
                            else
                                warning('does not meet a condition')
                            end
                        end
        

                    end

                    if strcmp(comparisonName, 'orientation_minus_baseline')
                        newMatrix(1, polarIndex, ri, si) = mean(horizontal);
                        newMatrix(2, polarIndex, ri, si) = mean(upperright);
                        newMatrix(3, polarIndex, ri, si) = mean(vertical);
                        newMatrix(4, polarIndex, ri, si) = mean(upperleft);
                    else
                        newMatrix(1, polarIndex, ri, si) = mean(rightwards);
                        newMatrix(2, polarIndex, ri, si) = mean(upperrightwards);
                        newMatrix(3, polarIndex, ri, si) = mean(upwards);
                        newMatrix(4, polarIndex, ri, si) = mean(upperleftwards);
                        newMatrix(5, polarIndex, ri, si) = mean(leftwards);
                        newMatrix(6, polarIndex, ri, si) = mean(lowerleftwards);
                        newMatrix(7, polarIndex, ri, si) = mean(downwards);
                        newMatrix(8, polarIndex, ri, si) = mean(lowerrightwards);
                    end
                
                   
                end
            end
        end

    else
        error('Must be project "da" to assume derivation of absolute directions.');
    end

end