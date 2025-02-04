function totalROIvertices = getROIidxs(subj, roiname, surfSize)
    lhSize = surfSize(1);
    lh_label = read_label(subj, fullfile('retinotopy_RE', sprintf('lh.%s_REmanual', roiname))); 
    rh_label = read_label(subj, fullfile('retinotopy_RE', sprintf('rh.%s_REmanual', roiname))); 
    totalROIvertices = [lh_label(:,1)+1; rh_label(:,1)+lhSize+1]; 
end