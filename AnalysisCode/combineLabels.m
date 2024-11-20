
cd('/Users/rje257/Documents/GitHub/DriftingGrating/AnalysisCode')

bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/derivatives/freesurfer';
setup_user('rania', bidsDir)

subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-0397', ...
    'sub-0442', 'sub-wlsubj121', ...
    'sub-wlsubj123', 'sub-wlsubj124', 'sub-wlsubj127', 'sub-0395', 'sub-0426', ...
    'sub-0427', 'sub-0250'}; 

for si=1:length(subjects)
    subjectname = subjects{si};
    labelDir = fullfile(bidsDir, subjectname, 'label', 'retinotopy_RE');
    cd(labelDir)

    [status, cmdout] = system(['mri_mergelabels -i lh.V2d_REmanual.label -i lh.V2v_REmanual.label -o lh.V2_REmanual.label'])

    [status, cmdout] = system(['mri_mergelabels -i rh.V2d_REmanual.label -i rh.V2v_REmanual.label -o rh.V2_REmanual.label'])


    [status, cmdout] = system(['mri_mergelabels -i lh.V3d_REmanual.label -i lh.V3v_REmanual.label -o lh.V3_REmanual.label'])

    [status, cmdout] = system(['mri_mergelabels -i rh.V3d_REmanual.label -i rh.V3v_REmanual.label -o rh.V3_REmanual.label'])
end