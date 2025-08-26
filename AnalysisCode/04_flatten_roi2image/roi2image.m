% plot label as image

%bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
bidsDir = '/Volumes/server/Projects/Project_dg/data_bids/';

addpath('/Users/rje257/Documents/GitHub/reToolbox/')


fslDir = '/Users/rje257/fsl/';
freesurferDir = '/Applications/freesurfer/7.4.1/';
githubDir = '/Users/rje257/Documents/GitHub/';

% FSL settings
PATH = getenv('PATH'); setenv('PATH', [PATH ':' fslDir '/bin']); % add freesurfer/bin to path
setenv('FSLDIR', fslDir);

% freesurfer settings
PATH = getenv('PATH'); setenv('PATH', [PATH ':' freesurferDir '/bin']); % add freesurfer/bin to path
setenv('FREESURFER_HOME', freesurferDir);
addpath(genpath(fullfile(freesurferDir, 'matlab')));
addpath(genpath(fullfile(freesurferDir, 'fsfast')));
setenv('SUBJECTS_DIR', [bidsDir '/derivatives/freesurfer']);

hemis = {'lh'}; %, 'rh'};
roinames = {'hMTcomplex'};
subjectlist =  {'sub-wlsubj127'};

for ss=1:numel(subjectlist)

    subject=subjectlist{ss};
    fsdir = fullfile(bidsDir, 'derivatives', 'freesurfer', subject);
    prfvistamovdir = retrieve_sesname(fullfile(bidsDir, 'derivatives', 'prfvista_mov', subject));

    for ri=1:numel(roinames)
    
        roiname = roinames{ri};
    
        roiIndx = [];
        for hi=1:numel(hemis)
    
            hh=hemis{hi};
    
            roi_label = fullfile(fsdir, 'label', 'retinotopy_RE', sprintf('%s.%s_REmanual.label', hh, roiname)); % _REmanual
            roiData = readmatrix(roi_label, 'FileType','text');
            
            roiIndx = [roiIndx; roiData(:,1)+1];

            prfvistamov_eccname = fullfile(prfvistamovdir, sprintf('%s.eccen.mgz', hh));
            eccVol = MRIread(prfvistamov_eccname); eccVol = eccVol.vol;
            eccVol = eccVol(roiIndx);

            % Extract RAS coordinates and label values
            x = roiData(:, 2); % X-coordinates in RAS
            y = roiData(:, 3); % Y-coordinates in RAS
            z = roiData(:, 4); % Z-coordinates in RAS (optional for 2D)
            labelValues = roiData(:, 5); % Label values

            scatter3(x, y, z, 50, eccVol, 'filled'); % 30 is the marker size
            colormap(jet); % Apply a colormap
            colorbar; % Display the colorbar
            caxis([0, 12]); % Set colorbar range from 0 to 12
            xlabel('X (RAS)');
            ylabel('Y (RAS)');
            zlabel('Z (RAS)');
            title('Label Plot (3D RAS)');

            % Create a 2D parameterized grid (e.g., flattening x-y plane)
            u = linspace(min(x), max(x), 100); % Define the grid in x
            v = linspace(min(y), max(y), 100); % Define the grid in y
            [U, V] = meshgrid(u, v);           % Create a 2D grid

            % Interpolate eccentricity values onto the 2D grid
            flattenedEccVol = griddata(x, y, eccVol, U, V, 'linear');


            % Display the flattened data
            figure;
            imagesc(u, v, flattenedEccVol); % Use u and v as axes
            colormap(jet);                  % Apply a colormap
            colorbar;                       % Show colorbar
            caxis([0, 12]); 
            xlabel('Flattened X');
            ylabel('Flattened Y');
            title('Flattened Surface with Eccentricity Values');
            
            % Adjust axis properties for visualization
            axis equal;
            axis tight;

        end
    end
end