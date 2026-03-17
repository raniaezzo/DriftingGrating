function [proConditions, conConditions, allConditions] = retrieveProConIdx(projectName, comparisonName, subset)

% recently edited to add vertical horizontal. Prior versions only included
% radialvstang for da, or "main cardinal" for da and dg

% these labels are the same because the betas were organized this way

% now subset=1 is vertical for DG and radial for DA (considered main /
% non-derived subset)

    % only used for conditions NOT derived (to index medianBOLDpa)
    if strcmp(comparisonName, 'motion_minus_orientation')
        if subset==1 && (strcmp(projectName, 'da') || (strcmp(projectName, 'dots')))
            proConditions = [9, 11]; 
            conConditions = [8, 10]; 
        elseif subset==1 && (strcmp(projectName, 'dg'))   % these are the same indices as da, based on how I organized the stimulus labels
            proConditions = [9, 11];                      % but I am making this explicit (checked)
            conConditions = [8, 10]; 
        else                        % same for dg, da
            proConditions = 8:11;
            conConditions = 12:15;
        end
        allConditions = [proConditions, conConditions];
    elseif strcmp(comparisonName, 'motion_minus_baseline')
        if subset==1 && (strcmp(projectName, 'da') || (strcmp(projectName, 'dots')))
            proConditions = 20:21;
            conConditions = 18:19;
        elseif subset==1 && (strcmp(projectName, 'dg'))
            proConditions = 20:21;                        % checked
            conConditions = 18:19;
        else                        % same for dg, da
            proConditions = 18:21;
            conConditions = 22:25;
        end
        allConditions = [proConditions, conConditions];
    elseif strcmp(comparisonName, 'orientation_minus_baseline')
        if subset==1 && (strcmp(projectName, 'da') || (strcmp(projectName, 'dots')))
            proConditions = 27;
            conConditions = 26;
        elseif subset==1 && (strcmp(projectName, 'dg'))
            proConditions = 27;                           % checked
            conConditions = 26;
        else                        % same for dg, da
            proConditions = 26:27;
            conConditions = 28:29;
        end
        allConditions = [proConditions, conConditions];
    end

end