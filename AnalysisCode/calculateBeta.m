function [y, yHat, betas, varexp] = calculateBeta(data, design, nEvents)

    % use canonical hRF
    hRF = getcanonicalhrf(1,1);

    % abs (to prevent flipping sign) - done in GLMestimatemodel
    mean_overT = mean(abs(data),2); 

    % doing PSC here - subtract 1 to remove baseline signal
    y = ((data./mean_overT)-1)*100;

    %y = data; % replaced with above to remove PSC
    
    eventMat = design(:,1:nEvents);

    pred_hRF = zeros(size(eventMat,1), 2);

     for i = 1:size(eventMat, 2)
        temp = conv(eventMat(:, i), hRF);
        % Truncate the convolution result to match the length of the original data
        pred_hRF(:, i) = temp(1:size(eventMat,1));
     end

    % Design: events, noise, constant, drift
    design = [pred_hRF, design(:,nEvents+1:end)];
    betas = pinv(design) * y';
    betas = betas';

    yHat = (design * betas')'; % this includes noise

    % calculate error (r^2)
    SSR = (((y - yHat).^2));
    SST = ((y - mean(y,2)).^2);
    varexp = 1 - (sum(SSR,2) ./ sum(SST,2));

end