function [dgTable, daTable] = reportROIvsHarmonicComparison(varargin)
% REPORTROIVSHARMONICCOMPARISON  Side-by-side ROI-binning vs harmonic-model
% asymmetry estimates for dgMatched7 and da (n=7 both), the same
% comparison summarized in Table S3, but printed for both projects with
% both methods' own point estimate, 68%/95% CI, and significance.
%
%   reportROIvsHarmonicComparison()
%   reportROIvsHarmonicComparison('roi', 'V2')
%
% ROI-binning estimates come from fitAsymmetryRegression.m's cache
% (regressionResults/<fitLabel>/<roi>.mat, mainCardinal/derivedCardinal/
% mainSubset/derivedSubset slot order); harmonic estimates come from
% fitHarmonicVertexModel.m's cache (regressionResultsHarmonic/<fitLabel>/
% <roi>.mat, project-invariant horizontal-vertical/cardinal-oblique/
% radial-tangential/polar-cardinal-oblique order). Both are remapped into
% the same project-invariant 4-asymmetry order here so the two methods'
% numbers for "radial minus tangential" (etc) sit in the same row for
% both dg and da, even though which raw regression slot that asymmetry
% occupies swaps between the two projects.

p = inputParser;
p.addParameter('roi', 'V1', @ischar);
p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/', @ischar);
p.parse(varargin{:});
opt = p.Results;

harmonicDisplayLabels = {'Horizontal minus Vertical','Cardinal minus Oblique','Radial minus Tangential','Polar Cardinal minus Polar Oblique'};

% bFromSlot(k) = which ROI-binning slot (mainCardinal=1/derivedCardinal=2/
% mainSubset=3/derivedSubset=4) holds harmonic asymmetry k (in
% harmonicDisplayLabels order) -- the inverse of fitHarmonicVertexModel.m's
% own slotFromB for each project.
bFromSlot.dg = [3, 1, 4, 2];
bFromSlot.da = [4, 2, 3, 1];

dgTable = oneProjectTable('dg', 'dgMatched7', bFromSlot.dg, harmonicDisplayLabels, opt);
daTable = oneProjectTable('da', 'da', bFromSlot.da, harmonicDisplayLabels, opt);

end

function T = oneProjectTable(projectName, fitLabel, bFromSlotThisProject, harmonicDisplayLabels, opt)

roiFile = fullfile(opt.bidsDir, 'derivatives', 'summaryTables', 'regressionResults', fitLabel, sprintf('%s.mat', opt.roi));
harmFile = fullfile(opt.bidsDir, 'derivatives', 'summaryTables', 'regressionResultsHarmonic', fitLabel, sprintf('%s.mat', opt.roi));
if ~isfile(roiFile)
    error('reportROIvsHarmonicComparison:missingROIFit', 'No fitAsymmetryRegression.m cache at %s -- run it first.', roiFile);
end
if ~isfile(harmFile)
    error('reportROIvsHarmonicComparison:missingHarmonicFit', 'No fitHarmonicVertexModel.m cache at %s -- run it first.', harmFile);
end
Froi = load(roiFile);
Fharm = load(harmFile);

roiEst = Froi.estimates(bFromSlotThisProject)';
roiCoeffs = Froi.coeffs(bFromSlotThisProject, :);
[roiCI68, roiCI95, roiSig] = ciAndSig(roiCoeffs);

harmEst = Fharm.harmonicEstimates';
harmCoeffs = Fharm.harmonicCoeffs;
[harmCI68, harmCI95, harmSig] = ciAndSig(harmCoeffs);

T = table(harmonicDisplayLabels', ...
    roiEst, roiCI68(:,1), roiCI68(:,2), roiCI95(:,1), roiCI95(:,2), roiSig, ...
    harmEst, harmCI68(:,1), harmCI68(:,2), harmCI95(:,1), harmCI95(:,2), harmSig, ...
    'VariableNames', {'Asymmetry', ...
    'ROI_Estimate', 'ROI_CI68_lower', 'ROI_CI68_upper', 'ROI_CI95_lower', 'ROI_CI95_upper', 'ROI_Sig', ...
    'Harmonic_Estimate', 'Harmonic_CI68_lower', 'Harmonic_CI68_upper', 'Harmonic_CI95_lower', 'Harmonic_CI95_upper', 'Harmonic_Sig'});

fprintf('%s (%s, n=%d): ROI-binning vs harmonic method / %s:\n', projectName, fitLabel, Froi.nSubj, opt.roi);
disp(T)

end

function [ci68, ci95, sig] = ciAndSig(bootDraws)
    ci68 = prctile(bootDraws, [16 84], 2);
    ci95 = prctile(bootDraws, [2.5 97.5], 2);
    sig = cell(size(bootDraws,1), 1);
    for hi = 1:size(bootDraws,1)
        if ci95(hi,1) > 0 || ci95(hi,2) < 0
            sig{hi} = '**';
        elseif ci68(hi,1) > 0 || ci68(hi,2) < 0
            sig{hi} = '*';
        else
            sig{hi} = '';
        end
    end
end
