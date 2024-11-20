function filteredPrfBins = retriveRetData(bidsDir, retFolder, subj, polarAngleBinWidth, minECC, maxECC, minVAREXP)

    polarAngles = 0:45:315;
    hemis = {'lh'; 'rh'};
    retDir = dir(fullfile(bidsDir, 'derivatives', retFolder, subj, '**/stimfiles.mat'));
    retDir = retDir.folder;

    for hi=1:numel(hemis)
            hemi = hemis{hi};
            ret_data.(sprintf('%s_pa', hemi)) = MRIread(fullfile(retDir, sprintf('%s.angle_adj.mgz', hemi)));
            ret_data.(sprintf('%s_ecc', hemi)) = MRIread(fullfile(retDir, sprintf('%s.eccen.mgz', hemi)));
            ret_data.(sprintf('%s_vexp', hemi)) = MRIread(fullfile(retDir, sprintf('%s.vexpl.mgz', hemi)));
            ret_data.(sprintf('%s_sigma', hemi)) = MRIread(fullfile(retDir, sprintf('%s.sigma.mgz', hemi)));      
    end

    % matrix: each column is vertex, rows are: pa, ecc, r^2, size
    prfData = [map_theta(ret_data.lh_pa.vol), map_theta(ret_data.rh_pa.vol); ...
        ret_data.lh_ecc.vol, ret_data.rh_ecc.vol ; ...
        ret_data.lh_vexp.vol, ret_data.rh_vexp.vol; ...
        ret_data.lh_sigma.vol, ret_data.rh_sigma.vol];

    filteredPrfBins = nan(1,length(prfData));

    % ultimately end with 8 numerical bins (0,45,90,135,180,225,315)
    % which is assigned to the wedge

    % and filter out by R^2 and eccentricity

    % loop through pa bins:
    for pa=1:numel(polarAngles)
        pangle = polarAngles(pa);
        
        % Calculate the lower and upper bounds of the angular bin
        lower_bound = mod(pangle - polarAngleBinWidth/2 + 180, 360) - 180;
        upper_bound = mod(pangle + polarAngleBinWidth/2 + 180, 360) - 180;

        if lower_bound<0
            lower_bound = 360+lower_bound;
        end
        if upper_bound<0
            upper_bound = 360+upper_bound;
        end
        
        % Find values within the specified angular bin
        if lower_bound>upper_bound
            currPA_bool = (prfData(1,:) >= lower_bound) | (prfData(1,:) <= upper_bound);
        else
            currPA_bool = (prfData(1,:) >= lower_bound) & (prfData(1,:) <= upper_bound);
        end
            
        % restrict eccentricities to be within a range
        currECC_bool = (prfData(2,:) >= minECC) & (prfData(2,:) <= maxECC);

        % restrict variance explained to be a certain minimum
        currVAREXP_bool = (prfData(3,:) >= minVAREXP);

        % if it meets the 3 criteria
        filterRet = currPA_bool+currECC_bool+currVAREXP_bool == 3;
        validIdx = find(filterRet>0 );

        filteredPrfBins(validIdx) = pangle; % assign a value (e.g., 0 degrees)

%         % find intersection between validIdx and label_idx
%         finalIdx = intersect(label_idx, validIdx);
    end
    

end