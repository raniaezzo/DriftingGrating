function convert_dgEDF(edfPath)

    datfilePath = strrep(edfPath, 'edf', 'dat');
    msgfilePath = strrep(edfPath, 'edf', 'msg');

    if ~isfile(datfilePath) || ~isfile(msgfilePath)

        % convert EDF to ASC
        sprintf('converting %s', edfPath)
    
        edf2ascFunction = '/usr/local/bin/edf2asc';
    
        [~,~] = system([edf2ascFunction,' ',edfPath,' -e -y']);
        movefile(strrep(edfPath, 'edf', 'asc'), msgfilePath); % rename part1 asc to msg (strs)
    
        [~,~] = system([edf2ascFunction,' ',edfPath,' -s -miss -1.0 -y']);
        movefile(strrep(edfPath, 'edf', 'asc'), datfilePath); % rename part2 asc to dat (#s)
    
        % do a check that the ASC is saved properly
        datContent = readmatrix(datfilePath);
    
        if size(datContent,2) > 4 % if there are more than 4 columns, parsing did not work
    
            numNonNaN = sum(~isnan(datContent), 2); % count how many rows do not have 4 valid values
            invalid_rows = find(numNonNaN ~= 4);
    
            % Process each invalid row
            for i = invalid_rows'
                % Set the entire row to NaN
                datContent(i, :) = NaN;
            
                % Ensure it's not the first row (avoid index error)
                if i > 1
                    datContent(i, 1) = datContent(i-1, 1) + 1; % Set column 1 as previous row +1
                else
                    datContent(i, 1) = 1; % If first row is invalid, set to 1
                end
            
                % Set columns 2 and 3 to -1, column 4 to 0
                datContent(i, 2:3) = -1;
                datContent(i, 4) = 0;
            end
    
            % Initialize a new matrix for cleaned data
            numRows = size(datContent, 1);
            clean_data = NaN(numRows, 4); % Preallocate for efficiency
    
            % Loop through rows and extract 4 non-NaN values
            for i = 1:numRows
                clean_data(i, :) = datContent(i, ~isnan(datContent(i, :))); % Extract non-NaN values
            end
    
            delete(datfilePath); % delete old dat file with unequal spacing
            writematrix(clean_data, datfilePath, 'Delimiter', 'tab');
    
        end
    end
    
end
