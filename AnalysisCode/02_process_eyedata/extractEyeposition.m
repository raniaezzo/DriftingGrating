function allDat = extractEyeposition(allDat, datfilePath, timestampTable, trialnumInit, trialsPerRun)

    currRun = floor(((trialnumInit+1)/trialsPerRun)+1);

    % Filter the table based on the run
    condition = timestampTable.runNum == currRun;
    filteredTable = timestampTable(condition, :);
    
    
    currDat = readmatrix(datfilePath);
    
    % in rare cases I need to downsample because collected at higher
    % sampling rate. Extract the rows corresponding to the first occurrence of each unique timestamp
    timestamps = currDat(:, 1);
    [~, first_indices] = unique(timestamps, 'first');
    currDat = currDat(first_indices, :);
    
    for ti=1:trialsPerRun
        % Find the first row where the first column equals startValue
        startIndex = find(currDat(:, 1) == filteredTable.sTrial(ti), 1, 'first');

        % Find the last row where the first column equals endValue
        endIndex = find(currDat(:, 1) == filteredTable.eTrial(ti), 1, 'last');

        if ~isempty(startIndex) && ~isempty(endIndex) && startIndex <= endIndex
            segmentTrialT = currDat(startIndex:endIndex, 1);
            segmentTrialX = currDat(startIndex:endIndex, 2);
            segmentTrialY = currDat(startIndex:endIndex, 3);
            segmentTrialP = currDat(startIndex:endIndex, 4);
            
            allDat(trialnumInit+ti, 1:length(segmentTrialX), 1) = segmentTrialT;
            allDat(trialnumInit+ti, 1:length(segmentTrialX), 2) = segmentTrialX;
            allDat(trialnumInit+ti, 1:length(segmentTrialY), 3) = segmentTrialY;
            allDat(trialnumInit+ti, 1:length(segmentTrialP), 4) = segmentTrialP;
        else
            warning('Fill trial with nans b/c not available in DAT')
        end
    end

end