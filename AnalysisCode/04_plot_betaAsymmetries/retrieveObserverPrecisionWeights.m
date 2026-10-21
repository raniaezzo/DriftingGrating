function w = retrieveObserverPrecisionWeights(subjects, roiname, precisionWeightsSource)
% RETRIEVEOBSERVERPRECISIONWEIGHTS  Per-observer precision weight for a
% given cortical area, mirroring retrieveObserverGainWeights.m's pattern
% but ROI-aware (precision/reliability genuinely varies by cortical area,
% unlike gain -- see computeObserverPrecisionWeights.m).
%
%   w = retrieveObserverPrecisionWeights(subjects, roiname, precisionWeightsSource)
%
% subjects               - cell array of subject IDs, in the order data is
%                           indexed by elsewhere (same convention as
%                           retrieveObserverGainWeights.m).
% roiname                - cortical area name (e.g. 'V1'), matching the
%                           'roi' column of precisionWeightsSource.
% precisionWeightsSource - [] (PLACEHOLDER: returns uniform weights of 1,
%                           i.e. a no-op, until the precision-weighting
%                           method is finalized), OR a table with columns
%                           'subject', 'roi', 'weight' (one row per
%                           (subject, cortical area) pair) once real
%                           precision weights are wired in.
%
% Returns a 1 x numel(subjects) row vector, in the same subject order.

if isempty(precisionWeightsSource)
    w = ones(1, numel(subjects));
    return
end

w = nan(1, numel(subjects));
for si = 1:numel(subjects)
    row = strcmp(precisionWeightsSource.subject, subjects{si}) & strcmp(precisionWeightsSource.roi, roiname);
    if ~any(row)
        error('retrieveObserverPrecisionWeights:missing', ...
            'No precision weight found for subject %s, roi %s', subjects{si}, roiname);
    end
    w(si) = precisionWeightsSource.weight(row);
end

end
