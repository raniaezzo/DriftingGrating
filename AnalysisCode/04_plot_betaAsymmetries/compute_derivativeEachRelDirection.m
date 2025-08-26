function newMatrix = compute_derivativeEachRelDirection(medianBOLDpa, projectSettings)

    % Check project name
    if strcmp(projectSettings.projectName, 'dg')
        
        projectName = projectSettings.projectName;
        comparisonName = projectSettings.comparisonName;
        contrasts_dict = projectSettings.contrasts_dict;
    
        %[proConditions, conConditions, allConditions] = retrieveProConIdx(projectName, comparisonName, radialvstang);
    
        % I think I always need to assume main cardinal, since radial / tang
        % for dg is defined ad hoc, outside of this function
        allConditions = retrieveRelDirIdx(projectName, comparisonName);
    
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
                        tangential = [];
                        radial = [];
                        cclock = [];
                        clock = [];
                    else
                        inwards = [];
                        outwards = [];
                        tangclock = [];
                        tangcclock = [];
                        outclock = [];
                        outcclock = [];
                        inclock = [];
                        incclock = [];
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
                
                        % label each local RF based on the relative dict (e.g. 90 = out; 270 = in etc.)
                        global_dir = motionValue;
                        localMotionDirs = deriveLocalMotionfromGlobal(global_dir, polarAngles);
                        currentlocalMotionDir = localMotionDirs(polarIndex);
    
                        % here orientations include 0-316 due to function
                        % conversion (deriveLocalMotionfromGlobal)
                        % NOTE: these values 0-316 are NOT in absolute
                        % motion terms -- they are RELATIVE (90 is out)
                        if strcmp(comparisonName, 'orientation_minus_baseline')
                            if ismember(currentlocalMotionDir, [0,180]) % locally tang (cc c)
                                tangential = [tangential, currentVal];
                            elseif ismember(currentlocalMotionDir, [90, 270]) % locally radial (out in)
                                radial = [radial, currentVal];
                            elseif ismember(currentlocalMotionDir, [135, 315])
                                cclock = [cclock, currentVal];
                            elseif ismember(currentlocalMotionDir, [45, 225])
                                clock = [clock, currentVal];
                            else
                                warning('does not meet a condition')
                            end
                        else
                            if ismember(currentlocalMotionDir, 270) % locally inwards
                                inwards = [inwards, currentVal];
                            elseif ismember(currentlocalMotionDir, 90) % locally outwards
                                outwards = [outwards, currentVal]; 
                            elseif ismember(currentlocalMotionDir, 0) % locally clockwise
                                tangclock = [tangclock, currentVal];
                            elseif ismember(currentlocalMotionDir, 180) % locally cclockwise
                                tangcclock = [tangcclock, currentVal]; 
                            elseif ismember(currentlocalMotionDir, 45)
                                outclock = [outclock, currentVal]; 
                            elseif ismember(currentlocalMotionDir, 135)
                                outcclock = [outcclock, currentVal];
                            elseif ismember(currentlocalMotionDir, 315)
                                inclock = [inclock, currentVal];
                            elseif ismember(currentlocalMotionDir, 225)
                                incclock = [incclock, currentVal];
                            else
                                warning('does not meet a condition')
                            end
                        end
        

                    end

                    if strcmp(comparisonName, 'orientation_minus_baseline')
                        newMatrix(1, polarIndex, ri, si) = mean(tangential);
                        newMatrix(2, polarIndex, ri, si) = mean(radial);
                        newMatrix(3, polarIndex, ri, si) = mean(cclock);
                        newMatrix(4, polarIndex, ri, si) = mean(clock);
                    else
                        newMatrix(1, polarIndex, ri, si) = mean(inwards);
                        newMatrix(2, polarIndex, ri, si) = mean(outwards);
                        newMatrix(3, polarIndex, ri, si) = mean(tangclock);
                        newMatrix(4, polarIndex, ri, si) = mean(tangcclock);
                        newMatrix(5, polarIndex, ri, si) = mean(outclock);
                        newMatrix(6, polarIndex, ri, si) = mean(outcclock);
                        newMatrix(7, polarIndex, ri, si) = mean(inclock);
                        newMatrix(8, polarIndex, ri, si) = mean(incclock);
                    end
                
                   
                end
            end
        end

    else
        error('Must be project "dg" to assume derivation of relative directions.');
    end

end