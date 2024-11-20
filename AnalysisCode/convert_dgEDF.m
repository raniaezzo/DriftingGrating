function convert_dgEDF(edfPath)

    % convert EDF to ASC
    sprintf('converting %s', edfPath)

    edf2ascFunction = '/usr/local/bin/edf2asc';

    [~,~] = system([edf2ascFunction,' ',edfPath,' -e -y']);
    msgfilePath = strrep(edfPath, 'edf', 'msg');
    movefile(strrep(edfPath, 'edf', 'asc'), msgfilePath); % rename part1 asc to msg (strs)

    [~,~] = system([edf2ascFunction,' ',edfPath,' -s -miss -1.0 -y']);
    datfilePath = strrep(edfPath, 'edf', 'dat');
    movefile(strrep(edfPath, 'edf', 'asc'), datfilePath); % rename part2 asc to dat (#s)
    
end
