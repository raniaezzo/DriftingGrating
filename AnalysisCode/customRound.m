function roundedValue = customRound(value)
    % Check if the value is positive or negative
    if value >= 0
        % Round up for positive values
        roundedValue = ceil(value * 1000) / 1000;
    else
        % Round down for negative values
        roundedValue = floor(value * 1000) / 1000;
    end
end