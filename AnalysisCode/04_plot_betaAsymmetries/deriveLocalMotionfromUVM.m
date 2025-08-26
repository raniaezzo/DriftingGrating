function localMotionDirs = deriveLocalMotionfromUVM(UVM_dir, polarAngles)

    % checked : 90; 270; 0; 180; 45; 135; 225; 315

    % which polarAngle contains the "cartesian motion" definition
    paWtrueDir = 90; % this does not change

    % unordered polar angles
    standardAngles = 0:45:315;

    % compute angular offset from UVM
    angularOffset = standardAngles-paWtrueDir;

    % direction at UVM
    localMotionDirs_tmp = repmat(UVM_dir, [1, length(polarAngles)]);

    localMotionDirs_tmp = localMotionDirs_tmp + angularOffset;

    % this just subtracts 360 if values to over 360
    localMotionDirs_tmp(localMotionDirs_tmp >= 360) = localMotionDirs_tmp(localMotionDirs_tmp >= 360) - 360;

    % if negative, add 360
    localMotionDirs_tmp(localMotionDirs_tmp < 0) = localMotionDirs_tmp(localMotionDirs_tmp < 0) + 360;

    % re-arrange according to the order of the provided polar angle array
    [~, idx] = ismember(standardAngles, polarAngles);
    localMotionDirs = localMotionDirs_tmp(idx);

end