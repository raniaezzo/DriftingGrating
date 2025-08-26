clear; clc; close all;

%% ===== Load video =====
%videoFile = '/Users/rje257/Desktop/movie_forward/forwardfixed_cycles.mp4';
%videoFile = '/Users/rje257/Desktop/movie_objectmotion/objectmotion.mp4';
videoFile = '/Users/rje257/Desktop/movie_eyemovements/eyemovements_contrastfix.mp4';
v = VideoReader(videoFile);

frames = [];
f = 1;
while hasFrame(v)
    rgb = readFrame(v);
    gray = rgb2gray(rgb);  % convert to grayscale
    frames(:,:,f) = im2double(gray); %#ok<SAGROW>
    f = f + 1;
end

[H,W,numFrames] = size(frames);

%% ===== Define central ROI =====
cropFraction = 0.5;  % <-- smaller than before (0.5 = 50% of width/height)

regionWidth  = round(cropFraction * W);
regionHeight = round(cropFraction * H);

% force it to be square by taking the smaller of the two
regionSize = min(regionWidth, regionHeight);
regionWidth = regionSize;
regionHeight = regionSize;

% center the square
xStart = round((W - regionWidth)/2);
yStart = round((H - regionHeight)/2);
xEnd = xStart + regionWidth - 1;
yEnd = yStart + regionHeight - 1;

% ===== Choose block size for ~12 samples across =====
numSamples = 8; %15; % adjust between 10-15 as desired
blockSize = round(regionWidth / numSamples);

% ===== Create grid centers =====
xCenters = (xStart + blockSize/2) : blockSize : (xEnd - blockSize/2);
yCenters = (yStart + blockSize/2) : blockSize : (yEnd - blockSize/2);
[Xg,Yg] = meshgrid(xCenters, yCenters);

%% ===== Motion estimation parameters =====
searchRange = 6; % max search in pixels

%% ===== Prepare figure and video writer =====
fig = figure('Position',[100 100 400 400]);
outputVideo = VideoWriter('output2.mp4','MPEG-4');
outputVideo.FrameRate = v.FrameRate;
open(outputVideo);

%% ===== Process frames =====
for t = 1:(numFrames-1)
    A = frames(:,:,t);
    B = frames(:,:,t+1);

    U = zeros(size(Yg));
    V = zeros(size(Yg));

    % ----- BLOCK MATCHING -----
    for by = 1:size(Yg,1)
        for bx = 1:size(Xg,2)
            cy = Yg(by,bx);
            cx = Xg(by,bx);
            y1 = round(cy - blockSize/2 + 1);
            y2 = round(cy + blockSize/2);
            x1 = round(cx - blockSize/2 + 1);
            x2 = round(cx + blockSize/2);
            if (y1<1 || x1<1 || y2>H || x2>W), continue; end

            patchA = A(y1:y2,x1:x2);
            if sum(patchA(:)) == 0
                continue; % skip empty regions
            end

            bestErr = inf; bestDx = 0; bestDy = 0;
            for dy = -searchRange:searchRange
                for dx = -searchRange:searchRange
                    ys = y1+dy; ye = y2+dy;
                    xs = x1+dx; xe = x2+dx;
                    if (ys<1 || xs<1 || ye>H || xe>W), continue; end
                    patchB = B(ys:ye,xs:xe);
                    err = sum((patchA(:)-patchB(:)).^2);
                    if err < bestErr
                        bestErr = err;
                        bestDx = dx;
                        bestDy = dy;
                    end
                end
            end

            % % normalize and scale for visibility
            % mag = sqrt(bestDx^2 + bestDy^2);
            % if mag > 0
            %     scaleFactor = 30; % make arrows longer for visibility
            %     U(by,bx) = (bestDx / mag) * scaleFactor;
            %     V(by,bx) = (bestDy / mag) * scaleFactor;
            % end

            scaleFactor = 10;
            U(by,bx) = bestDx * scaleFactor;  % displacement in pixels between frames
            V(by,bx) = bestDy * scaleFactor;

        end
    end

    U = imgaussfilt(U,1); % smooth across neighboring blocks
    V = imgaussfilt(V,1);

    N = 15; % average over 15 frames
    idx = mod(t-1, N) + 1;
    U_hist(:,:,idx) = U;
    V_hist(:,:,idx) = V;

    % ----- PLOT RESULTS -----
    % subplot(1,2,1);
    % imshow(B(yStart:yEnd, xStart:xEnd), []);
    % %title(sprintf('Frame %d (ROI)', t+1));
    % axis square;
    % 
    % % subplot(1,2,2);
    % % quiver(Xg,Yg,U,V,0,'Color','b');
    % % axis ij; axis([xStart xEnd yStart yEnd]); axis square;
    % % title('Local motion directions');
    % 
    % subplot(1,2,2);
    set(gca, 'Color', 'w');    % axes background white
    set(gcf, 'Color', 'w');    % figure background white
    set(gca, 'XColor', 'w', 'YColor', 'w');  % make axis lines and ticks white
    box on;                                  % ensure the box is visible
    if t >= N % average over 15 frames
        U_avg = mean(U_hist, 3);
        V_avg = mean(V_hist, 3);
    

        quiver(Xg, Yg, U_avg, V_avg, 0, 'k', 'LineWidth', 2);
        axis ij;
        axis([xStart-20 xEnd-20 yStart-20 yEnd-20]);
        axis square;
        xticks([]); yticks([]); 
        set(gca, 'Color', 'w');    % axes background white
        set(gcf, 'Color', 'w');    % figure background white
        set(gca, 'XColor', 'w', 'YColor', 'w');  % make axis lines and ticks white
        box on; 
        %title(sprintf('Average motion over last %d frames', N));
    end

    drawnow;

    % ----- Save frame to output video -----
    frame = getframe(fig);
    writeVideo(outputVideo, frame);
end

close(outputVideo);
disp('Movie saved as output.avi');