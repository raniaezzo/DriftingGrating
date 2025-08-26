function writeMGZfile(bidsDir, subjectname, ses, data, derivativesFolder, fileName)

    derivativesSubfolder = fullfile(derivativesFolder, subjectname, ses);
    
    if ~isfolder(derivativesSubfolder)
        mkdir(derivativesSubfolder)
    end

    hh = {'lh','rh'};
    % get dimensions
    hSize = get_surfsize(subjectname); % get hemi size
    hSizeIdx=[1,hSize(1);hSize(1)+1,sum(hSize)];
    % save mgz
    origfile = fullfile(bidsDir, 'derivatives','freesurfer', subjectname, 'mri', 'orig.mgz');
    mgz = MRIread(origfile);
    
    for hi=1:numel(hh)
        mgz.vol = data(hSizeIdx(hi,1):hSizeIdx(hi,2), 1);
        savePath = fullfile(derivativesSubfolder, sprintf('%s.%s.mgz',hh{hi},fileName));
        MRIwrite(mgz, savePath);
    end

end