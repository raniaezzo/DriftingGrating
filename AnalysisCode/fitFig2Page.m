function gcf_edit = fitFig2Page(gcf_edit)

% Get the figure's current position and aspect ratio
    pos = get(gcf_edit, 'Position'); % [left, bottom, width, height]
    aspectRatio = pos(3) / pos(4); % width/height
    
    % Define landscape page size (in inches, default for PDF)
    pageWidth = 11; % Width of US letter in landscape
    pageHeight = 8.5; % Height of US letter in landscape
    
    % Adjust the figure size while maintaining aspect ratio
    if aspectRatio > pageWidth / pageHeight
        % Wider than the page: scale based on width
        scaledWidth = pageWidth;
        scaledHeight = pageWidth / aspectRatio;
    else
        % Taller than the page: scale based on height
        scaledHeight = pageHeight;
        scaledWidth = pageHeight * aspectRatio;
    end
    
    % Set the paper size and position
    set(gcf_edit, 'PaperUnits', 'inches');
    set(gcf_edit, 'PaperSize', [pageWidth, pageHeight]);
    set(gcf_edit, 'PaperPosition', [(pageWidth - scaledWidth) / 2, ... % Center horizontally
                               (pageHeight - scaledHeight) / 2, ... % Center vertically
                               scaledWidth, scaledHeight]); % Scaled dimensions

    end