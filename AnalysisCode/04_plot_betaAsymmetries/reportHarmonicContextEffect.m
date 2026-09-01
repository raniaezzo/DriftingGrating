function contextTable = reportHarmonicContextEffect(varargin)
% REPORTHARMONICCONTEXTEFFECT  Paired dg-vs-da context effect for each of
% the 4 harmonic asymmetries, from fitHarmonicVertexModel.m's cached
% dgMatched7/da fits. Prints a table and returns it.
%
%   reportHarmonicContextEffect()
%   reportHarmonicContextEffect('roi', 'V2')
%
% Uses dgMatched7 (NOT dg_all) because it is the same 7 subjects, in the
% same order, as da -- both fitHarmonicVertexModel.m calls build subjects
% from the identical literal list, so perObserverEstimates row i is the
% same subject in both caches, for every ROI. That is what makes the
% bootstrap below genuinely PAIRED: each draw resamples one set of
% subject indices and applies it to BOTH projects' perObserverEstimates
% in the same iteration (same convention as
% fitAsymmetryRegression_dgVsDa.m's paired bootstrap), rather than
% subtracting two independently-bootstrapped distributions, which would
% inflate the apparent uncertainty with between-subject noise common to
% both experiments.
%
% perObserverEstimates (and therefore the context effect here) is in the
% project-invariant harmonicTermNames order (horizontal-vertical,
% cardinal-oblique, radial-tangential, polar-cardinal-oblique) -- the
% same 4 physical concepts in both projects, unlike the
% mainCardinal/derivedCardinal/mainSubset/derivedSubset slots whose
% concept meaning swaps between dg and da.

p = inputParser;
p.addParameter('roi', 'V1', @ischar);
p.addParameter('nBoot', 1000, @isnumeric);
p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/', @ischar);
p.parse(varargin{:});
opt = p.Results;

harmonicDisplayLabels = {'Horizontal minus Vertical','Cardinal minus Oblique','Radial minus Tangential','Polar Cardinal minus Polar Oblique'};

dgFile = fullfile(opt.bidsDir, 'derivatives', 'summaryTables', 'regressionResultsHarmonic', 'dgMatched7', sprintf('%s.mat', opt.roi));
daFile = fullfile(opt.bidsDir, 'derivatives', 'summaryTables', 'regressionResultsHarmonic', 'da', sprintf('%s.mat', opt.roi));
if ~isfile(dgFile) || ~isfile(daFile)
    error('reportHarmonicContextEffect:missingFit', ...
        'Cached harmonic fit not found for dgMatched7 and/or da / %s -- run fitHarmonicVertexModel(''dg'',''dgSubjectMode'',''matched'') and fitHarmonicVertexModel(''da'') first.', opt.roi);
end
Fdg = load(dgFile);
Fda = load(daFile);

if ~isequal(Fdg.subjects, Fda.subjects)
    error('reportHarmonicContextEffect:subjectMismatch', ...
        'dgMatched7 and da subject lists/order differ for %s -- the paired bootstrap below assumes row-for-row correspondence and cannot proceed.', opt.roi);
end
nSubj = Fdg.nSubj;

dgEst = Fdg.perObserverEstimates; % nSubj x 4, project-invariant order
daEst = Fda.perObserverEstimates;

% Each project's OWN asymmetry CI: read straight from
% fitHarmonicVertexModel.m's own cached bootstrap (harmonicCoeffs, 4 x
% nBoot, same project-invariant order) -- not recomputed here, so this
% matches exactly what that function printed for each project on its own.
[dgCI68, dgCI95, dgSig] = harmonicCIandSig(Fdg.harmonicCoeffs);
[daCI68, daCI95, daSig] = harmonicCIandSig(Fda.harmonicCoeffs);

pointEstimate = mean(dgEst, 1, 'omitnan') - mean(daEst, 1, 'omitnan'); % dg minus da -- context effect

rng(1, 'twister');
diffBoot = nan(4, opt.nBoot);
validObsIdx = find(~any(isnan(dgEst),2) & ~any(isnan(daEst),2));
for b = 1:opt.nBoot
    bIdx = validObsIdx(randi(numel(validObsIdx), numel(validObsIdx), 1)); % SAME resampled subjects applied to both projects, same draw
    diffBoot(:,b) = mean(dgEst(bIdx,:), 1)' - mean(daEst(bIdx,:), 1)';
end
[ci68, ci95, sig] = harmonicCIandSig(diffBoot);

contextTable = table(harmonicDisplayLabels', ...
    mean(dgEst,1,'omitnan')', dgCI68(:,1), dgCI68(:,2), dgCI95(:,1), dgCI95(:,2), dgSig, ...
    mean(daEst,1,'omitnan')', daCI68(:,1), daCI68(:,2), daCI95(:,1), daCI95(:,2), daSig, ...
    pointEstimate', ci68(:,1), ci68(:,2), ci95(:,1), ci95(:,2), sig, ...
    'VariableNames', {'Asymmetry', ...
    'dgMatched7', 'dg_CI68_lower', 'dg_CI68_upper', 'dg_CI95_lower', 'dg_CI95_upper', 'dg_Sig', ...
    'da', 'da_CI68_lower', 'da_CI68_upper', 'da_CI95_lower', 'da_CI95_upper', 'da_Sig', ...
    'ContextEffect_dgMinusDa', 'Context_CI68_lower', 'Context_CI68_upper', 'Context_CI95_lower', 'Context_CI95_upper', 'Context_Sig'});

fprintf('Paired dg-vs-da context effect (%d matched subjects) / %s:\n', nSubj, opt.roi);
disp(contextTable)

end

function [ci68, ci95, sig] = harmonicCIandSig(bootDraws)
% HARMONICCIANDSIG  68%/95% percentile CI and significance stars from a
% 4 x nBoot bootstrap draw matrix -- same convention used throughout this
% pipeline (** = 95% CI excludes 0, * = only 68% does).
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
