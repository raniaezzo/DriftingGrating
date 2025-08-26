function [proConditions, conConditions, allConditions] = retrieveProConIdx(projectName, comparisonName, radialvstang)

    % only used for conditions NOT derived (to index medianBOLDpa)
    if strcmp(comparisonName, 'motion_minus_orientation')
        if radialvstang==1 && strcmp(projectName, 'da')
            proConditions = [9, 11]; 
            conConditions = [8, 10]; 
        else                        % same for dg, da
            proConditions = 8:11;
            conConditions = 12:15;
        end
        allConditions = [proConditions, conConditions];
    elseif strcmp(comparisonName, 'motion_minus_baseline')
        if radialvstang==1 && strcmp(projectName, 'da')
            proConditions = 20:21;
            conConditions = 18:19;
        else                        % same for dg, da
            proConditions = 18:21;
            conConditions = 22:25;
        end
        allConditions = [proConditions, conConditions];
    elseif strcmp(comparisonName, 'orientation_minus_baseline')
        if radialvstang==1 && strcmp(projectName, 'da')
            proConditions = 27;
            conConditions = 26;
        else                        % same for dg, da
            proConditions = 26:27;
            conConditions = 28:29;
        end
        allConditions = [proConditions, conConditions];
    end

end