function [subjectyMin, subjectyMax] = getSubjectTTaxis(subj, plotValue)

    if strcmp(subj, 'sub-0395')
        subjectyMin = -1;
        subjectyMax = 2;
        scaleFactor = 0.6;
    elseif strcmp(subj, 'sub-0250')
        subjectyMin = -0.75;
        subjectyMax = 1.5;
        scaleFactor = 0.65;
    elseif strcmp(subj, 'sub-0255')
        subjectyMin = -0.5;
        subjectyMax = 1.5;
        scaleFactor = 0.6; 
    elseif strcmp(subj, 'sub-0201')
        subjectyMin = -0.5;
        subjectyMax = 0.6;
        scaleFactor = 0.7;
    elseif strcmp(subj, 'sub-0427')
        subjectyMin = -1.5;
        subjectyMax = 2.5;
        scaleFactor = 0.7;
    elseif strcmp(subj, 'sub-0037')
        subjectyMin = -0.6;
        subjectyMax = 0.8;
        scaleFactor = 0.85;
    elseif strcmp(subj, 'sub-0397')
        subjectyMin = -0.4;
        subjectyMax = 1;
        scaleFactor = 0.7;
    elseif strcmp(subj, 'sub-0426')
        subjectyMin = -0.5;
        subjectyMax = 1.5;
        scaleFactor = 0.7;
    elseif strcmp(subj, 'sub-0442')
        subjectyMin = -1;
        subjectyMax = 2;
        scaleFactor = 0.7;
    elseif strcmp(subj, 'sub-wlsubj121')
        subjectyMin = -0.6;
        subjectyMax = 1;
        scaleFactor = 0.75;
    elseif strcmp(subj, 'sub-wlsubj127')
        subjectyMin = -1;
        subjectyMax = 2;
        scaleFactor = 0.5;
    elseif strcmp(subj, 'sub-wlsubj123')
        subjectyMin = -0.75;
        subjectyMax = 1.75;
        scaleFactor = 0.7;
    elseif strcmp(subj, 'sub-wlsubj124')
        subjectyMin = -0.75;
        subjectyMax = 1.75;
        scaleFactor = 0.8;
    elseif strcmp(subj, 'allsubjects')
        subjectyMin = -0.5;
        subjectyMax = 1.25;  
        scaleFactor = 0.7;
    else
        subjectyMin = -1;
        subjectyMax = 1;
    end

    if strcmp(plotValue, 'subtractOri')
        subjectyMin = subjectyMin*scaleFactor;
        subjectyMax = subjectyMax*scaleFactor;
    end

end