function allDat_dva = removeShortsegs(allDat_dva, minDur)
    % Loop through each row (trial)
    for row = 1:size(allDat_dva, 1)
        % Get the data for the current row
        dataRow = allDat_dva(row, :);

        % Find indices where data is **not** NaN
        validIdx = find(~isnan(dataRow));

        if isempty(validIdx)
            continue; % Skip if the entire row is NaN
        end

        % Identify segment start and end
        segmentStart = validIdx([true, diff(validIdx) > 1]); % First element or discontinuity
        segmentEnd = validIdx([diff(validIdx) > 1, true]);   % Last element before discontinuity

        % Loop through each segment
        for seg = 1:length(segmentStart)
            segmentLength = segmentEnd(seg) - segmentStart(seg) + 1;
            if segmentLength < minDur
                % Set short segments to NaN
                dataRow(segmentStart(seg):segmentEnd(seg)) = NaN;
            end
        end

        % Store back into the main matrix
        allDat_dva(row, :) = dataRow;
    end
end