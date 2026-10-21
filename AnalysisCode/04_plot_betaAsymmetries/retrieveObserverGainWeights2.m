function gainWeights = retrieveObserverGainWeights2(subjects, roiname, gainWeightsSource)
% RETRIEVEOBSERVERGAINWEIGHTS2  Per-observer, per-cortical-area pRF gain,
% mirroring retrieveObserverPrecisionWeights.m's pattern (ROI-aware
% lookup from a long-format table). Supersedes retrieveObserverGainWeights.m,
% which returned one V1-only value per observer regardless of ROI.
%
%   gainWeights = retrieveObserverGainWeights2(subjects, roiname, gainWeightsSource)
%
% subjects        - cell array of subject IDs, in the order data is
%                    indexed by elsewhere (same convention as
%                    retrieveObserverGainWeights.m /
%                    retrieveObserverPrecisionWeights.m).
% roiname          - cortical area name (e.g. 'V1'), matching the 'roi'
%                    column of gainWeightsSource.
% gainWeightsSource - table with columns 'subject', 'roi', 'weight' (one
%                    row per (subject, cortical area) pair), e.g.
%                    gainTable from computeObserverGainWeightsByROI.m
%                    (loaded from gainSummaryByROI.mat). Unlike
%                    retrieveObserverPrecisionWeights.m's precisionWeightsSource,
%                    there is no [] no-op placeholder here: gain
%                    correction is always applied (it rescales each
%                    subject's raw data, not an optional weighting), so a
%                    missing source is an error rather than a silent 1s
%                    fallback.
%
% Returns a 1 x numel(subjects) row vector, in the same subject order. A
% subject with no row for this ROI (e.g. pMT/pMST, where a couple of
% subjects' retinotopy label/vertices don't survive the pRF inclusion
% criteria -- see computeObserverGainWeightsByROI.m's log) gets NaN, NOT
% an error: NaN propagates through subjectScale and the resulting
% gain-corrected value, and gets dropped as a missing row by fitlm
% downstream, same as any other subject with no usable data in a given
% ROI. Callers computing groupGain = exp(mean(log(gainWeights))) MUST use
% 'omitnan', or one missing subject would silently NaN out every other
% subject's correction factor for that ROI.
%
% See also COMPUTEOBSERVERGAINWEIGHTSBYROI, RETRIEVEOBSERVERPRECISIONWEIGHTS

if isempty(gainWeightsSource)
    error('retrieveObserverGainWeights2:missingSource', ...
        ['gainWeightsSource is empty -- gain correction is not optional (unlike precision, ' ...
         'it always rescales the data). Run computeObserverGainWeightsByROI.m and pass its ' ...
         'gainTable (e.g. via projectSettings.gainWeightsSource).']);
end

gainWeights = nan(1, numel(subjects));
for si = 1:numel(subjects)
    row = strcmp(gainWeightsSource.subject, subjects{si}) & strcmp(gainWeightsSource.roi, roiname);
    if ~any(row)
        continue % gainWeights(si) stays NaN -- see note above
    end
    gainWeights(si) = gainWeightsSource.weight(row);
end

end
