function localMotionDirs = deriveLocalMotionfromGlobal(global_dir, polarAngles)

    % this keeps the meaning of 90 (outward) but defined locally so:
    % 180: cc relative to gaze ; 0: clock relative to gaze ect.

    % manually checked : 90, 0, 45, 135, 180, 225, 270, 315

    % unordered polar angles
    standardAngles = 0:45:315;

    % because the added angle goes from UVM in the clockwise direction,
    % I create this array to re-order this starting PA 0-325
    % this is computed because the futher away from the UVM the more
    % angular offset. Instead of +- I use larger values that 180
    directionRefOffset = [90, 45, 0, 315, 270, 225, 180, 135];

    % direction at UVM
    localMotionDirs_tmp = repmat(global_dir, [1, length(polarAngles)]);

    localMotionDirs_tmp = localMotionDirs_tmp + directionRefOffset;

    % this just subtracts 360 if values to over 360
    localMotionDirs_tmp(localMotionDirs_tmp >= 360) = localMotionDirs_tmp(localMotionDirs_tmp >= 360) - 360;

    % if negative, add 360
    localMotionDirs_tmp(localMotionDirs_tmp < 0) = localMotionDirs_tmp(localMotionDirs_tmp < 0) + 360;

    % re-arrange according to the order of the provided polar angle array
    [~, idx] = ismember(standardAngles, polarAngles);
    localMotionDirs = localMotionDirs_tmp(idx);

end