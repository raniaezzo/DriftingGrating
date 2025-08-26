function finalTable = findMSs(allDatdva, finalTable, msLabels, metaDataTable)
    
    for trial=1:height(finalTable)
        runN = finalTable.runNum(trial);
        samplingRate = metaDataTable.SamplingRate(find(metaDataTable.RunID==runN));
    
        % go through each trial in allDat to extract indices with
        % microsaccades:
    
        % 1) filter x and y for detection
        idxs = allDatdva(trial,:,1); % EDF time
        x = allDatdva(trial,:, 2); % 
        y = allDatdva(trial,:, 3); % 

        % this fill fail in cases where the array ends with nans
        % and the nans were used as filler to ensure than all trials
        % are EXACTLY the same length. So now we will remove those nans
        % and optionally add them back later.
        [x, n_xendNans] = removeTrailingNaNs(x);
        [y, n_yendNans] = removeTrailingNaNs(y);

        sprintf('trial: %i',trial)

        if trial==105
            disp('stop')
        end

        if ~isempty(x) && ~isempty(y) % skip if the whole trial is invalid
            x_fil=filtfilt(fir1(35,0.05),1,x);
            y_fil=filtfilt(fir1(35,0.05),1,y);
        
            % position vector
            d=[x_fil,y_fil];
        
            % velocity vector
            v = computevelocity(d,samplingRate); 
        
            % return output cell array
            % first cell contains non-filtered saccades: no amplitude thresh;
                % and includes amplitudes < 0.01 and > 1 degree
            % second cell contains filtered microsaccades: amplitude thresh
            output = detectSacc(d, v, samplingRate);
    
            for ll=1:length(msLabels)
                if ~isempty(output{1})
                    if strcmp(msLabels{ll}, 'sMSacc')
                        timeIDs = idxs(output{1}(:,1)');
                    elseif strcmp(msLabels{ll}, 'eMSacc')
                        timeIDs = idxs(output{1}(:,2)');
                    end
                    finalTable.(msLabels{ll})(trial) = {timeIDs};
                end
            end
        end
        
    end

end

%%
function [trimmedArray, numTrailingNaNs] = removeTrailingNaNs(arr)
    % Ensure the input is a column vector
    arr = arr(:); 
    
    % Find the last non-NaN index
    lastValidIndex = find(~isnan(arr), 1, 'last');

    % Check if the array ends in NaNs
    if isempty(lastValidIndex) % If the entire array is NaN
        numTrailingNaNs = length(arr);
        trimmedArray = [];
    else
        numTrailingNaNs = length(arr) - lastValidIndex; % Count trailing NaNs
        trimmedArray = arr(1:lastValidIndex); % Remove them
    end
end
