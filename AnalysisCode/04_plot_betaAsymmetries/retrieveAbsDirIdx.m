function [allConditions] = retrieveAbsDirIdx(projectName, comparisonName)

    % contrastnames = {
    % 1- cardMsep
    % 2- oblMsep
    % 3 - allmValls
    % 4 - allsVblank
    % 5 - allmVblank
    % 6 - cardmVblank
    % 7 - oblmVblank
    % 8 - m0_v_s90
    % 9 - m90_v_s0
    % 10 - m180_v_s90
    % 11 - m270_v_s0
    % 12 - m45_v_s135
    % 13 - m135_v_s45
    % 14 - m225_v_s135
    % 15 - m315_v_s45
    % 16 - cardsVblank
    % 17 - oblsVblank
    % 18 - m0_v_b
    % 19 - m180_v_b 
    % 20 - m90_v_b 
    % 21 - m270_v_b 
    % 22 - m45_v_b 
    % 23 - m225_v_b 
    % 24 - m135_v_b
    % 25 - m315_v_b 
    % 26 - s0_v_b 
    % 27 - s90_v_b 
    % 28 - s45_v_b 
    % 29 - s135_v_b

    % (to index medianBOLDpa)
    % for derived values (e.g., DA), further computations are done later
    if strcmp(comparisonName, 'motion_minus_orientation')
        rightwards = [8]; 
        upperrightwards = [12]; 
        upwards = [9];
        upperleftwards = [13];
        leftwards = [10];
        lowerleftwards = [14];
        downwards = [11];
        lowerrightwards = [15];
        allConditions = [rightwards, upperrightwards, upwards, upperleftwards, ...
            leftwards, lowerleftwards, downwards, lowerrightwards];
    elseif strcmp(comparisonName, 'motion_minus_baseline')
        rightwards = [18]; 
        upperrightwards = [22]; 
        upwards = [20];
        upperleftwards = [24];
        leftwards = [19];
        lowerleftwards = [23];
        downwards = [21];
        lowerrightwards = [25];
        allConditions = [rightwards, upperrightwards, upwards, upperleftwards, ...
            leftwards, lowerleftwards, downwards, lowerrightwards];
    elseif strcmp(comparisonName, 'orientation_minus_baseline')
        horizontal = [26];
        upperright = [28];
        vertical = [27];
        upperleft = [29];
        allConditions = [horizontal, upperright, vertical, upperleft];
    end

end