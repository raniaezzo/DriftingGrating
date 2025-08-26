function [allConditions] = retrieveRelDirIdx(projectName, comparisonName)

    % contrastnames = {
    % 1- cardMsep (in, out, clock, cclock)
    % 2- oblMsep (other - not in, out, clock, cclock)
    % 3 - allmValls
    % 4 - allsVblank
    % 5 - allmVblank
    % 6 - cardmVblank (in, out, clock, cclock)
    % 7 - oblmVblank (other - not in, out, clock, cclock)
    % 8 - m0_v_s90 (clock)
    % 9 - m90_v_s0 (outward)
    % 10 - m180_v_s90 (cclock)
    % 11 - m270_v_s0 (inward)
    % 12 - m45_v_s135 (out clock)
    % 13 - m135_v_s45 (out cclock)
    % 14 - m225_v_s135 (in cclock)
    % 15 - m315_v_s45 (in clock)
    % 16 - cardsVblank (pinwheel, annuli)
    % 17 - oblsVblank (spirals)
    % 18 - m0_v_b (clock)
    % 19 - m180_v_b (cclock)
    % 20 - m90_v_b (out)
    % 21 - m270_v_b (in)
    % 22 - m45_v_b (out clock)
    % 23 - m225_v_b (in cclock)
    % 24 - m135_v_b (out cclock)
    % 25 - m315_v_b (in clock)
    % 26 - s0_v_b (annulus)
    % 27 - s90_v_b (pinwheel)
    % 28 - s45_v_b (clock)
    % 29 - s135_v_b (cclock)

    % (to index medianBOLDpa)
    % for derived values (e.g., DG), further computations are done later
    if strcmp(comparisonName, 'motion_minus_orientation')
        inwards = [11]; 
        outwards = [9]; 
        tangclock = [8];
        tangcclock = [10];
        outclock = [12];
        outcclock = [13];
        inclock = [15];
        incclock = [14];
        allConditions = [inwards, outwards, tangclock, tangcclock, ...
            outclock, outcclock, inclock, incclock];
    elseif strcmp(comparisonName, 'motion_minus_baseline')
        inwards = [21]; 
        outwards = [20]; 
        tangclock = [18];
        tangcclock = [19];
        outclock = [22];
        outcclock = [24];
        inclock = [25];
        incclock = [23];
        allConditions = [inwards, outwards, tangclock, tangcclock, ...
            outclock, outcclock, inclock, incclock];
    elseif strcmp(comparisonName, 'orientation_minus_baseline')
        tangential = [26];
        radial = [27];
        cclock = [29];
        clock = [28];
        allConditions = [tangential, radial, cclock, clock];
    end

end