function [predicted, R2, betas] = refitGroupMeanTrace(meanObserved, designMatrix, stimdur_s, tr_s)
% refitGroupMeanTrace - re-fit a fresh GLM to a cross-subject group-mean
% ROI-average time series (can't reuse individual-subject GLMsingle
% per-vertex HRF fits once traces are averaged across subjects).
%
% <designMatrix> is the run's condition-onset matrix (nTRs x nConditions,
% one column per condition, 1 at onset TR). Each column is convolved with
% the canonical HRF; an intercept column is added; betas are estimated by
% ordinary least squares.
%
% Returns the predicted time series, R^2 (1 - SSres/SStot) between
% <meanObserved> and <predicted>, and the fitted betas (last one is the
% intercept).

nTRs = size(designMatrix, 1);
nCond = size(designMatrix, 2);

hrf = getcanonicalhrf(stimdur_s, tr_s);
hrf = hrf(:);

X = zeros(nTRs, nCond);
for ci = 1:nCond
    convolved = conv(designMatrix(:, ci), hrf);
    X(:, ci) = convolved(1:nTRs);
end
X = [X, ones(nTRs, 1)]; % intercept

meanObserved = meanObserved(:);
betas = X \ meanObserved;
predicted = X * betas;

ssRes = sum((meanObserved - predicted).^2);
ssTot = sum((meanObserved - mean(meanObserved)).^2);
R2 = (1 - ssRes / ssTot) * 100;

predicted = predicted';

end
