function dataTable = extractTimestamps(msgfilePath, dataTable, trialnumInit, trialsPerRun)

    currRun = floor(((trialnumInit+1)/trialsPerRun)+1);

    msgfid = fopen(msgfilePath,'r');

    % Read the file contents into a character array
    fileData = fread(msgfid, '*char')';
    fclose(msgfid);

    % Split the data into lines
    % Convert the character array to a string (if using newer versions of MATLAB, use "string" function)
    fileData = string(fileData);
    lines = strsplit(fileData, '\n'); % Split the data into lines

    % Initialize arrays to store the data

    currTrial = trialnumInit;

    % Loop through each line to extract the relevant timestamps and labels
    for i = 1:length(lines)
        line = strtrim(lines{i});

        % trial start 
        if contains(line, 'TRIAL') && contains(line, 'START') && currTrial==trialnumInit
            currTrial = currTrial+1;
            tokens = regexp(line, 'MSG\s+(\d+)\s+TRIAL (\d+) START', 'tokens');
            if ~isempty(tokens)
                dataTable.sTrial(currTrial) = str2double(tokens{1}{1});
                dataTable.sStim(currTrial) = str2double(tokens{1}{1});
            end

        % stimulus offset
        elseif contains(line, 'STIMULUS OFFSET')  && currTrial>trialnumInit
            tokens = regexp(line, 'MSG\s+(\d+)\s+STIMULUS OFFSET', 'tokens');
            if ~isempty(tokens)
                dataTable.eStim(currTrial) = str2double(tokens{1}{1});
            end

        % start / end of saccade
         elseif startsWith(line, 'ESACC')  && currTrial>trialnumInit
            tokens = regexp(line, '^ESACC\s+\w+\s+(\d+)\s+(\d+)', 'tokens');
            if ~isempty(tokens)
                sTimestamps = str2double(tokens{1}{1});
                eTimestamps = str2double(tokens{1}{2});
                dataTable.sSacc{currTrial} = [dataTable.sSacc{currTrial}, sTimestamps];
                dataTable.eSacc{currTrial} = [dataTable.eSacc{currTrial}, eTimestamps];
            end

         % start / end of blink
         elseif startsWith(line, 'EBLINK')  && currTrial>trialnumInit
            tokens = regexp(line, '^EBLINK\s+\w+\s+(\d+)\s+(\d+)', 'tokens');
            if ~isempty(tokens)
                sTimestamps = str2double(tokens{1}{1});
                eTimestamps = str2double(tokens{1}{2});
                dataTable.sBlink{currTrial} = [dataTable.sBlink{currTrial}, sTimestamps];
                dataTable.eBlink{currTrial} = [dataTable.eBlink{currTrial}, eTimestamps];
            end         

         % end of trial
         elseif contains(line, 'TRIAL') && contains(line, 'END')  && currTrial>trialnumInit
            tokens = regexp(line, 'MSG\s+(\d+)\s+TRIAL (\d+) END', 'tokens');
            if ~isempty(tokens)
                dataTable.eTrial(currTrial) = str2double(tokens{1}{1});
                dataTable.trialNum(currTrial) = currTrial-trialnumInit;
                dataTable.runNum(currTrial) = currRun;
            end

            if currTrial==trialnumInit+trialsPerRun
                disp('ending')
                break;
            else
                currTrial = currTrial+1;
                dataTable.sTrial(currTrial) = str2double(tokens{1}{1});
                dataTable.sStim(currTrial) = str2double(tokens{1}{1});
            end

        end
    end

    % if eStim == NaN , it is a blank trial
    blankTrials = isnan(dataTable.eStim); % Find the rows where column 'eStim' is NaN
    % Set the corresponding values in column 'B' to NaN
    dataTable.sStim(blankTrials) = NaN;

end
