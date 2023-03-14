function mask = createmask(outdoorimage, radius)

    % compute dimensions of image (height, width, and colors)
    [imageH, imageW, imageC] = size(outdoorimage);
    
    % initialize mask matrix of all zeros
    mask = ones(imageH, imageW);

    [y x] = size(mask); %define y,x as size of array
    r = radius; %imageH/2; %define radius of a circle
    for i=1:y
        for j=1:x
            if ((i-y/2)^2)+((j-x/2)^2)<(r^2)  %define origin is at the center
                mask(i,j) = nan;  %define array inside the circle eq. = 1
            end
        end
    end 
end