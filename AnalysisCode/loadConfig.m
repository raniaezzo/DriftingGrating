function projectSettings = loadConfig(githubDir)

    roisfile = fullfile(githubDir, 'DriftingGrating', 'AnalysisCode', 'jsons', 'ROIS.json');
    str = fileread(roisfile); rois_dict = jsondecode(str);
    rois = {rois_dict.ROIs.filename};
    
    % read in file for axes limits
    axesLimfile = fullfile(githubDir, 'DriftingGrating', 'AnalysisCode', 'jsons', 'polaraxes_limits.json');
    str = fileread(axesLimfile);
    axes_limits = jsondecode(str);
        % EXAMPLE USE
        % ax.RLim = [axes_limits.(projectName).(comparisonName).ROIs_early.min ...
        %             axes_limits.(projectName).(comparisonName).ROIs_early.max]

    % read in file for axes limits
    axesLimfile = fullfile(githubDir, 'DriftingGrating', 'AnalysisCode', 'jsons', 'pairwise_limits.json');
    str = fileread(axesLimfile);
    pairaxes_limits = jsondecode(str);

    % read in file for axes limits for equal weight pairwise plots
    axesLimfile = fullfile(githubDir, 'DriftingGrating', 'AnalysisCode', 'jsons', 'pairwise_equalPAweight_limits.json');
    str = fileread(axesLimfile);
    pairaxes_PAew_limits = jsondecode(str);

    % read in file for axes limits
    ttaSubjectsLimfile = fullfile(githubDir, 'DriftingGrating', 'AnalysisCode', 'jsons', 'subjectwise_tta_limits.json');
    str = fileread(ttaSubjectsLimfile);
    subjectwise_tta_limits = jsondecode(str);
    
    colorsfile = fullfile(githubDir, 'DriftingGrating', 'AnalysisCode', 'jsons', 'COLORS.json');
    str = fileread(colorsfile);
    colors_data = jsondecode(str);
        % EXAMPLE USE
        % color_values_pro = colors_data.conditions.(projectName).("mainCardinalVsMainOblique").color_pro;
        % color_values_con = colors_data.conditions.(projectName).("mainCardinalVsMainOblique").color_con;
    
    contrastsfile = fullfile(githubDir, 'DriftingGrating', 'AnalysisCode', 'jsons', 'CONTRASTS.json');
    str = fileread(contrastsfile); contrasts_dict = jsondecode(str);
        % EXAMPLE USE
        % contrastnames = {contrasts_dict.contrasts.(strcat(projectName, '_contrast_name'))};

    projectSettings.gitDir = githubDir;
    projectSettings.rois = rois;
    projectSettings.axes_limits = axes_limits;
    projectSettings.pairaxes_limits = pairaxes_limits;
    projectSettings.pairaxes_PAew_limits = pairaxes_PAew_limits;
    projectSettings.colors_data = colors_data;
    projectSettings.contrasts_dict = contrasts_dict;
    projectSettings.subjectwise_tta_limits = subjectwise_tta_limits;

end

