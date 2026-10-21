function excludedRunIdx = findDuplicateDesignRuns(matrices_onset, nRunsCheck)
% findDuplicateDesignRuns - detect runs whose condition-onset design
% matrix is identical to an earlier run's for the same subject (i.e. the
% same trial randomization was used/shown more than once), which are
% therefore not statistically independent repeats.
%
% Found in practice for da/sub-wlsubj124: S02_design_Run1.mat and
% S02_design_Run2.mat (experimentalOutput/da/02/) have byte-identical
% trialMat content -- confirmed not a coincidence (52 trials, multiple
% randomized conditions) and already baked into matrices_onset as loaded
% from rawInfo.mat. A cross-run data-vs-model correlation check showed the
% real BOLD data for both runs matches this shared design well, i.e. the
% same stimulus sequence really was presented twice (not a mislabeling
% bug) -- but averaging both into a group mean would double-count that one
% trial sequence relative to every other run.
%
% For each pair of runs (i, j) with i < j <= nRunsCheck whose
% matrices_onset{i} and matrices_onset{j} are identical (NaN-safe), the
% LATER run (j) is marked excluded, keeping the earlier occurrence.
%
% Returns excludedRunIdx, a row vector of run indices to exclude (may be
% empty).

nRunsUse = min(numel(matrices_onset), nRunsCheck);
excludedRunIdx = [];
for i = 1:nRunsUse
    if ismember(i, excludedRunIdx)
        continue % already excluded as someone else's duplicate
    end
    for j = (i+1):nRunsUse
        if ismember(j, excludedRunIdx)
            continue
        end
        if isequaln(matrices_onset{i}, matrices_onset{j})
            excludedRunIdx(end+1) = j; %#ok<AGROW>
        end
    end
end

end
