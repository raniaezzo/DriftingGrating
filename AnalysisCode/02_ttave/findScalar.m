function fitted_advMotionVals = findScalar(grandMean, advMotionVals)

%         % Compute the scaling factor k
%         k = (advMotionVals * grandMean') / (grandMean * grandMean');
%         % Scale array2
%         fitted_advMotionVals = k * advMotionVals;

        % Add a column of ones for the intercept
        X = [ones(length(grandMean), 1), grandMean'];
        beta = regress(advMotionVals', X);

        % Extract the intercept and scaling factor
        intercept = beta(1); scaling_factor = beta(2);
        % Scale array2 using GLM
        fitted_advMotionVals = intercept + scaling_factor * grandMean;
end
