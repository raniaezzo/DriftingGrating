

numRhoLines = 8;
numThetaLines = 16;
maxRadius = 1;

figure('Color', 'w'); hold on; axis equal;
thetaVals = linspace(0, 2*pi, 1000);

% Draw concentric circles (rho lines)
for r = linspace(maxRadius/numRhoLines, maxRadius, numRhoLines)
    x = r * cos(thetaVals);
    y = r * sin(thetaVals);
    plot(x, y, 'k-', 'LineWidth',2);
end

% % Draw radial lines (theta lines)
% thetaLines = linspace(0, 2*pi, numThetaLines+1);
% for t = thetaLines(1:end-1)
%     plot([0 cos(t)] * maxRadius, [0 sin(t)] * maxRadius, 'k-', 'LineWidth',2);
% end
set(gca, 'Visible', 'off');  
axis([-maxRadius maxRadius -maxRadius maxRadius]);

figure('Color', 'w'); hold on; axis equal;
numXLines = 8;
numYLines = 8;
maxValue = 1;

xVals = linspace(-maxValue, maxValue, numXLines + 1);
yVals = linspace(-maxValue, maxValue, numYLines + 1);

% Draw vertical grid lines
for x = xVals
    plot([x x], [-maxValue, maxValue], 'k-', 'LineWidth',2);
end

% Draw horizontal grid lines
for y = yVals
    plot([-maxValue, maxValue], [y y], 'k-', 'LineWidth',2);
end

axis([-maxValue maxValue -maxValue maxValue]);
set(gca, 'Visible', 'off');  

%%

numSpirals = 10; 
numPoints = 1000; 
maxRadius = 1; 
spiralTightness = .5;  % Higher = more spiraling

figure('Color', 'w');
hold on; axis equal off;

startAngles = linspace(0, 2*pi, numSpirals+1);
startAngles(end) = []; % remove duplicate endpoint

for startTheta = startAngles
    theta = linspace(0, spiralTightness * pi, numPoints);
    rho = linspace(0, maxRadius, numPoints);

    % Offset the spiral by the start angle
    x = rho .* cos(theta + startTheta);
    y = rho .* sin(theta + startTheta);

    plot(x, y, 'k-', 'LineWidth',2);
end

axis([-maxRadius maxRadius -maxRadius maxRadius]);