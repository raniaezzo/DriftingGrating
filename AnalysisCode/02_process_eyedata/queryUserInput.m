function userInputs = queryUserInput(xPos, yPos, trialTimeStamps)
% Loops through n and queries user input while displaying plots.
%
% INPUT:
%   xPos - An (n x m) matrix where each row is plotted sequentially.
%
% OUTPUT:
%   userInputs - A cell array (1 x n), where each cell contains:
%       - A user-defined vector (if provided)
%       - An empty array if only ENTER was pressed

    n = size(xPos, 1); % Number of iterations
    userInputs = cell(n,1);
    %userInputsStart = cell(n, 1); % Preallocate cell array
    %userInputsEnd = cell(n, 1); % Preallocate cell array

    for tt = 1:n
        % Plot the current row of xPos_Hmotion
        figure(1); clf;
        set(gcf, 'Position', [18 846 1778 491])
        plot(xPos(tt, :), 'b-'); % Plot current row
        hold on
        plot(yPos(tt, :), 'r-'); % Plot current row
        hold on
        yline(0, '--', 'Color', [.5 .5 .5])
        legend({'x', 'y'})
        title(sprintf('Iteration %d of %d', tt, n));
        xlabel('Time since trial onset (ms)'); ylabel('distance from center (deg)');
        ylim([-2, 2])
        % Set major ticks every 1,000 and minor ticks every 100
        xlim([1, size(xPos,2)]) % Ensure the range is correct
        xticks(0:1000:size(xPos,2)) % Major ticks at every 1,000
        xticklabels(xticks) % Keep labels only for major ticks
        ax = gca;
        ax.XAxis.MinorTickValues = 0:100:size(xPos,2); % Minor ticks every 100
        ax.XMinorGrid = 'on'; % Optional: Show minor grid
        drawnow;

        figure(2); clf;
        set(gcf, 'Position', [15 275 620 491])
        xlim([-3 3])
        ylim([-3, 3])
        fixationSize = getMarkerSize(0.36);
        scatter(0,0,fixationSize, 'or', 'filled', 'LineWidth', 10)
        hold on
        plot(xPos(tt,:), yPos(tt,:), 'b-', 'LineWidth', 1); % 'b-' specifies a blue line
        hold on
        % Overlay semi-transparent scatter points
        scatter(xPos(tt, :), yPos(tt, :), 5, 'k', 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeAlpha', 0.01);
        ylabel('deg', 'FontSize', 20)
        xlabel('deg', 'FontSize', 20)
        xlim([-3 3])
        ylim([-3, 3])
        ax = gca; % Get the current axis
        ax.FontSize = 20; 
        drawnow;

        edfIdxArray = trialTimeStamps(tt, :);

        % Prompt user for input
        i = 1; % Initialize index for storage
        userInputsStart = [];
        userInputsEnd = [];
        while true
            % Display input prompt
            fprintf('Enter a START time (1 to 3000), or press ENTER to exit:\n');
        
            % Set up key press detection
            set(gcf, 'KeyPressFcn', @adjustLimits); % Set keypress function
        
            % Wait for user input while allowing zooming
            userInputStart = input('', 's');  % Console remains active for input
        
            % Check if user pressed ENTER to exit
            if isempty(userInputStart)
                break; % Exit the loop
            end
        
            % Convert input to numeric
            parsedStart = str2double(userInputStart) + 1;
        
            % Validate input
            if isnan(parsedStart) || parsedStart < 0 || parsedStart > length(xPos) || mod(parsedStart, 1) ~= 0
                disp('Invalid input. Please enter a whole number between 0 and length(xPos).');
                continue; % Restart the loop
            end
        
            % Ask for END input
            while true
                userInputEnd = input(sprintf('Enter an END index (> %d and ≤ 3000): ', parsedStart), 's');
                parsedEnd = str2double(userInputEnd) + 1;
        
                % Validate END input
                if isnan(parsedEnd) || parsedEnd <= parsedStart || parsedEnd > length(xPos) || mod(parsedEnd, 1) ~= 0
                    disp('Invalid input. END must be greater than START and within valid range.');
                else
                    break; % Valid END input, exit loop
                end
            end
        
            % Store values
            userInputsStart = [userInputsStart, edfIdxArray(parsedStart)];
            userInputsEnd = [userInputsEnd, edfIdxArray(parsedEnd)];
            i = i + 1; % Increment storage index
        end
        

    userInputs{tt, 1} = userInputsStart;
    userInputs{tt, 2} = userInputsEnd;
    disp('User input collection complete.');
    
    
    hold off;

    end

end

function markerSize = getMarkerSize(diameter)
    ax = gca; % Get current axis
    
    % Get axis limits
    xl = ax.XLim;
    yl = ax.YLim;
    
    % Get figure pixel size of the axis
    axUnits = get(ax, 'Units');
    set(ax, 'Units', 'pixels'); % Temporarily change to pixels
    axPos = get(ax, 'Position'); % Get axis size in pixels
    set(ax, 'Units', axUnits); % Restore original units

    % Compute data units per pixel
    xScale = diff(xl) / axPos(3); % X-axis scale (data units per pixel)
    yScale = diff(yl) / axPos(4); % Y-axis scale (data units per pixel)
    scale = mean([xScale, yScale]); % Average to preserve aspect ratio

    % Convert diameter to scatter marker area in points²
    markerSize = (diameter / scale) ^ 2; % Scatter marker size is area in pt²
end

% Function to adjust zooming without breaking user input
function adjustLimits(~, event)
    % Get axes for figure(1) and figure(2)
    fig1 = figure(1); ax1 = gca;
    fig2 = figure(2); ax2 = gca;

    fig2Children = get(ax2, 'Children');
    fixationHandle = fig2Children(3); 

    % Handle UP arrow (double limits)
    if strcmp(event.Key, 'uparrow')
        set(0, 'CurrentFigure', fig1);
        %xlim(ax1, xlim(ax1) * 2);
        ylim(ax1, ylim(ax1) * 2);

        set(0, 'CurrentFigure', fig2);
        xlim(ax2, xlim(ax2) * 2);
        ylim(ax2, ylim(ax2) * 2);

        % Scale fixation size down
        fixationSize = get(fixationHandle, 'SizeData');
        set(fixationHandle, 'SizeData', fixationSize / 2);

    % Handle DOWN arrow (halve limits)
    elseif strcmp(event.Key, 'downarrow')
        set(0, 'CurrentFigure', fig1);
        %xlim(ax1, xlim(ax1) / 2);
        ylim(ax1, ylim(ax1) / 2);

        set(0, 'CurrentFigure', fig2);
        xlim(ax2, xlim(ax2) / 2);
        ylim(ax2, ylim(ax2) / 2);

        % Scale fixation size up
        fixationSize = get(fixationHandle, 'SizeData');
        set(fixationHandle, 'SizeData', fixationSize * 2);
    end
end
