% will convert design to include only stim ONSETs (required for
% GLMsingle/GLMdenoise)
function cellArray = filterOnset(cellArray)
    for i = 1:numel(cellArray)
        % Get the matrix from the current cell unit
        matrix = cellArray{i};
        
        % Initialize a new matrix with the same size, filled with zeros
        new_matrix = zeros(size(matrix));
        
        % Loop through each column (condition)
        for j = 1:size(matrix, 2)
            % Find all 1s in the current column
            event_indices = find(matrix(:, j) == 1);
            
            % Loop through the event indices to find the onset of each event
            for k = 1:length(event_indices)
                if k == 1 || event_indices(k) > event_indices(k-1) + 1
                    % Set the onset (first 1 of an event) in the new matrix
                    new_matrix(event_indices(k), j) = 1;
                end
            end
        end
        
        % Overwrite the original matrix with the modified matrix
        cellArray{i} = new_matrix;
    end
end