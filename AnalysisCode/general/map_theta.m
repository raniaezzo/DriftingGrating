function theta = map_theta(theta)
    % function to convert Benson deg to conventional deg (note: these values are VISUAL FIELD)
    
    theta = theta*-1;
    
    orig_size = size(theta);
    theta = reshape(theta, 1, []);

    for tt=1:length(theta)
        
        if theta(tt) >= -90 && theta(tt) < 0
        % Shift the range (-90, 0) to (0, 90)
        	theta(tt) = theta(tt) + 90;
            
        elseif theta(tt) >= 0 && theta(tt) < 90
        % Shift the range (0, 90) to (90, 180)
            theta(tt) = theta(tt) + 90;

        elseif theta(tt) >= 90 && theta(tt) <= 180
        % Shift the range (90, 180) to (180, 270)
            theta(tt) = theta(tt) + 90;
            
        elseif theta(tt) >= -180 && theta(tt) < -90
        % Shift the range (-180, -90) to (270, 360)
            theta(tt) = theta(tt) + 450;
        end
        
    end
    
    theta = reshape(theta, orig_size);
end