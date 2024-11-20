function [designfiles,dg_id, dg_numID, run] = format_desmats(bidsDir, designDir, sub, TR_s)

    designfiles = {};

    % retrieve alternate ID
    participant_file = struct2table(tdfread(fullfile(bidsDir,'participants.tsv')));
    rowInd = contains(cellstr(participant_file.participant_id), sub);
    dg_id = participant_file.alternate_id(rowInd,:);
    dg_numID = regexp(dg_id, '\d+', 'match'); dg_numID = dg_numID{1};

    disp('DG UNIQUE ID:')
    dg_numID

    % compute how many runs exist for this subject
    logFile = fullfile(designDir, dg_numID, 'runlog.txt');
    fid = fopen(logFile, 'r');
    runs = fscanf(fid, '%s');
    runarray = split(runs, 'Run');
    totalRuns = str2num(runarray{end});
    run = 1:totalRuns;

    % set order of design matrix
    % (1) M_0, 
    % (2) M_90, 
    % (3) M_180, 
    % (4) M_270, 
    % (5) M_45, 
    % (6) M_135, 
    % (7) M_225, 
    % (8) M_315, 
    % (9) S_0, 
    % (10) S_90, 
    % (11) S_45, 
    % (12) S_135, 
    % (13) B
    
    % condition_labels = {'M_card_d0_o90','M_card_d90_o0','M_card_d180_o90','M_card_d270_o0', ...
    %     'M_obl_d45_o135','M_obl_d135_o45','M_obl_d225_o135','M_obl_d315_o45', ...
    %     'S_card_o0', 'S_card_o90', 'S_obl_o45', 'S_obl_o135', 'B'};
    condition_labels = {'M_0','M_90','M_180','M_270', ...
        'M_45','M_135','M_225','M_315', ...
        'S_0', 'S_90', 'S_45', 'S_135', 'B'};
    condition_dmPos = 1:13;
    mapCond = containers.Map(condition_labels, condition_dmPos);

    for ri=run % run is the array 1:X

        % convert matrix from experimental script for GLM
    
        path_to_dmat = fullfile(designDir, dg_numID, sprintf('Run%i', ri), ...
            sprintf('%s_design_Run%i.mat', dg_id, ri));
        load(path_to_dmat, 'expDes')
        
        n_trs = expDes.runPadding_s+ ...
            (expDes.stimDur_s+expDes.itiDur_s)*expDes.nb_trials+ ...
            expDes.runPadding_s;
        
        % unique conditions (regressors)
        n_conds = expDes.nb_trials/expDes.nb_repeat;
        
        design = zeros(n_trs, n_conds);
        
        % padding period (in TRs)
        runPadding_trs = expDes.runPadding_s/TR_s;
        % stimulus period (in TRs)
        stimDur_trs = expDes.stimDur_s/TR_s;
        % iti period (in TRs)
        itiDur_trs = expDes.itiDur_s/TR_s;
        
        %%
        currentTR=runPadding_trs+1;
        for i=1:expDes.nb_trials
        
            stimMat = zeros(stimDur_trs,n_conds);
            isiMat = zeros(itiDur_trs,n_conds);
        
            shortV = expDes.trialMat(i,:);
            
            if shortV(2) == 0 % baseline
                level = 'B';
                dir = '';
            elseif shortV(2) == 1 % static
                level = 'S_';
                dir = num2str(shortV(3));
            elseif shortV(2) == 2 % motion
                level = 'M_';
                dir = num2str(shortV(4));
            end
        
            condition = strcat(level,dir);
            idx = mapCond(condition);
        
            stimMat(:,idx) = 1;
        
            trialMat = [stimMat; isiMat];
        
            endTR = currentTR+4;
        
            design(currentTR:endTR,:) = trialMat;
        
            currentTR = endTR+1;
           
        end
        
        designfiles{ri} = design;

    end

end