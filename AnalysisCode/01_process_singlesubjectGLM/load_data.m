function datafiles = load_data(bidsDir,task,space,fileType,sub,ses,runs)

% Inputs:
%     bidsDir:        Path to project folder containing the derivative folder -e.g. '/Volumes/Vision/MRI/Decoding'
%     task:           task name from the scans -e.g. 'Cue'
%     space:          data space -e.g. 'fsnative', 'fsaverage6', etc
%     fileType:       data file type -e.g. 'gii', 'mgh', 'nii', etc
%     sub:            Subject ID -e.g. 'sub-0201'
%     ses:            Session ID -e.g. 'ses-01'
%     runs:            Run ID -e.g. [1:10]

% Outputs:
%     datafiles:        data stored in 1 by #runs cell matrix

if isempty(ses)
    subDir = sprintf('%s/derivatives/fmriprep/%s/func',bidsDir,sub);
else
    subDir = sprintf('%s/derivatives/fmriprep/%s/%s/func',bidsDir,sub,ses);
end


switch space % switch between volumn or surface
    
    case 'T1w'
        
        datafiles = cell(1,length(runs));
        
        for iRun = 1:length(runs)
            
            whichRun = runs(iRun);
            
            disp(['Loading: ' sprintf('%s/*%s*%s_*%s_desc-preproc_bold%s',subDir,task,num2str(whichRun),space,fileType)]);
            tmpDir = dir(sprintf('%s/*task-%s_run-%s_*%s_desc-preproc_bold%s*',subDir,task,num2str(whichRun),space,fileType));
            datafiles{iRun} = niftiread([subDir '/' tmpDir.name]);
            
        end
        
        datafiles(cellfun(@isempty,datafiles)) = [];
        
    case {'fsnative','fsaverage','fsaverage5','fsaverage6'}
        
        hemi = {'L';'R'};
        datafiles = cell(1,length(runs)); % initialize for all the runs
        
        for iRun = 1:length(runs)
            
            whichRun = runs(iRun); % current run #
            func = cell(2,1); % initialize for 2 hemi
            
            for iH = 1:numel(hemi)
                
                %fileName = sprintf('%s/derivatives/fmriprep/%s/%s/func/%s_%s_task-%s_run-%s_space-fsnative_hemi-%s_bold.func',bidsDir,sub,ses,sub,ses,task,num2str(whichRun),hemi{iH});
                %input = [fileName '.gii'];

                inputStrcut = dir(sprintf('%s/derivatives/fmriprep/%s/%s/func/%s_%s_task-%s*_run-%s_space-fsnative_hemi-%s_bold.func.gii',bidsDir,sub,ses,sub,ses,task,num2str(whichRun),hemi{iH}));
                input = fullfile(inputStrcut.folder, inputStrcut.name);

                %output = [fileName fileType]; % the file type that we want to load

                output = strrep(input, '.gii', fileType);
                
                % check to see if data exists in the desired fileType, if not,
                % mir_convert file from gii
                if ~exist(output, 'file') && exist(input, 'file')
                    disp(['File does not exist in ' fileType ' format, converting from .gii ...'])
                    system(['mri_convert ' input ' ' output]);
                elseif ~exist(output, 'file') && ~exist(input, 'file')
                    error(['File does not exist in neither .gii or ' fileType ' format. Check name.'])
                end
                
                disp(['Loading: ' output])
                
                switch fileType
                    case '.gii'
                        tmp  = gifti(output);
                        func{iH}  = func.cdata;
                    case {'.mgh','.mgz'}
                        tmp = MRIread(output);
                        func{iH} = squeeze(tmp.vol);
                    otherwise
                        error('file type not valid')
                end
                datafiles{iRun} = cat(1,func{:});
                
                
            end
        end
        
    otherwise
        error('data space not valid')
end

end

%%