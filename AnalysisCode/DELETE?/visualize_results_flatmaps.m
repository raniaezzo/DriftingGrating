% clear workspace
clearvars;
close all;
clc;

% ses-01 (dg)
% ses-02 (retinotopy)
% ses-03 (localizer)

%% Setup A
% check if in correct working directory
[ParentFolderPath, Parentfolder] = fileparts(pwd);
if ~strcmp(Parentfolder, 'retinotopy_visualization')
    disp('WARNING: You are not in the correct directory. Go into retinotopy_visualization.')
end

% setup variables & paths 
[projectDir, githubDir, freesurferDir] = setup_user_visualization('rania');
subject           = '0201';
dataStr           = 'fsnative'; % fsaverage; fsaverage6; fsnative
modelType         = 'fine';
projectName       = 'dg';
localizerName     = 'loc';
retinotopyName    = 'Retinotopy_NYUAD';

addpath(genpath(pwd));

ses               = 'ses-02'; % for retinotopy 
prfOptsPathCoarse = fullfile(projectDir, 'derivatives', 'stim_apertures', sprintf('sub-%s', subject), sprintf('%s', ses), 'prfOptsCoarse.json');
[averageScans,stimwidthdeg,opt] = getPRFOpts(prfOptsPathCoarse);
stimradiusdeg = stimwidthdeg*0.5;

if ~strcmp(dataStr, 'fsnative')
    sub = dataStr;
else
    sub = strcat('sub-',subject);
end


%% Setup B: Exact brain atlases -> Wang and Glasser in fsaverage and fsnative space:

% check that Glasser & Wang atlases exist. If not create all the labels
% this is already inside the github repo -> https://github.com/WinawerLab/atlasmgz

returndir = pwd;

% Generate wang and glasser atlases (in all directories)
glasser_dir = fullfile(projectDir, 'derivatives', 'freesurfer', sub, 'label', 'Glasser2016');
if ~isfolder(glasser_dir)
    cd(fullfile(githubDir, 'atlasmgz'))
    system(['sh createAtlasLabels.sh ' subject ' ' projectDir])
end

cd(returndir);

%% %% Setup C: Define Localizer Design and Analysis Params

ses               = 'ses-03'; % for retinotopy 
run = 1:8;
rois = {'mt','mstL','mstR','fst'}; % later omit fst
stimdur_s = 10;
tr_s = 1;
hrf_opt = 'assume';
resampling = 0; % means fit fully (don't bootstrap or cross-validate)


%% PART I:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FUNCTIONAL LOCALIZER%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Description:
% This runs a basic GLM where each unique event (moving-static) is a unique regressor
% (1) moving, 
% (2) static, 
% Produces 2 beta maps (one for each regressor).

% Step 1) create design matrix
tmp{1} = repmat(repelem([0 1;1 0],10,1),15,1);
matrices = repmat(tmp(1), 1, max(run));

%% load in the functional localizer data

% load_data (function in wpToolbox) - convert to .mgh if not already
datafiles = load_data(projectDir,localizerName,dataStr,'.mgh',sub,ses,run);

% run glm (from GLMdenoise repo)
for iRoi = 1:numel(rois) % leave out FST
    results.(char(rois{iRoi})) = GLMestimatemodel(matrices([iRoi iRoi+4]),datafiles([iRoi iRoi+4]), ...
        stimdur_s,tr_s,hrf_opt,[],resampling);
end

%% Contrast maps

% Produces contrast maps (voxelwise):
% Estimated amplitude for motion - estimated amplitude for static

betamap_mt = results.mt.modelmd{1,2};
betamap_mstL = results.mstL.modelmd{1,2};
betamap_mstR = results.mstR.modelmd{1,2};

% contrast map for moving - static
contrasts.mt_motionVstatic = betamap_mt(:,1) - betamap_mt(:,2);
contrasts.mtL_motionVstatic = betamap_mstR(:,1) - betamap_mstR(:,2);
contrasts.mtR_motionVstatic = betamap_mstL(:,1) - betamap_mstL(:,2);

% convert to correct space (work around until all run with same recon-all)
% first retrieve left / right dimensions
vals = betamap_mt(:,1); % entire brain dimensions
setenv('SUBJECTS_DIR', [projectDir '/derivatives/freesurfer']); 
hSize = get_surfsize(sub); % get hemi size
hSizeIdx=[1,hSize(1);hSize(1)+1,numel(vals)];  

% for MT center
lhcontrast = double(contrasts.mt_motionVstatic(hSizeIdx(1,1):hSizeIdx(1,2)));  
rhcontrast = double(contrasts.mt_motionVstatic(hSizeIdx(2,1):hSizeIdx(2,2)));  
% convert to fsaverage space
% fsaverage_lhcontrast = cvntransfertosubject(sub,'fsaverage', lhcontrast, 'lh', 'nearest', 'orig', 'orig');
% fsaverage_rhcontrast = cvntransfertosubject(sub,'fsaverage', rhcontrast, 'rh', 'nearest', 'orig', 'orig');

% fsaverageresults.contrasts.mtcenter_motionVstatic = [fsaverage_lhcontrast; fsaverage_rhcontrast];
results.contrasts.mtcenter_motionVstatic = [lhcontrast; rhcontrast];


% take the motion from left stimulus for right hemi (& visa versa)
% then combine them:
% for MT periphery
lhcontrast = double(contrasts.mtL_motionVstatic(hSizeIdx(1,1):hSizeIdx(1,2)));  
rhcontrast = double(contrasts.mtR_motionVstatic(hSizeIdx(2,1):hSizeIdx(2,2)));  
% convert to fsaverage space
% fsaverage_lhcontrast = cvntransfertosubject(sub,'fsaverage', lhcontrast, 'lh', 'nearest', 'orig', 'orig');
% fsaverage_rhcontrast = cvntransfertosubject(sub,'fsaverage', rhcontrast, 'rh', 'nearest', 'orig', 'orig');
% 
% fsaverageresults.contrasts.mtperiphery_motionVstatic = [fsaverage_lhcontrast; fsaverage_rhcontrast];
% 
% fsaverageresults.contrasts.mtall_motionVstatic = fsaverageresults.contrasts.mtcenter_motionVstatic + ...
%     fsaverageresults.contrasts.mtperiphery_motionVstatic;

results.contrasts.mtperiphery_motionVstatic = [lhcontrast; rhcontrast];

results.contrasts.mtall_motionVstatic = results.contrasts.mtcenter_motionVstatic + ...
    results.contrasts.mtperiphery_motionVstatic;

%results.contrasts = contrasts;

%visualize R2 for each roi
%close all;
% mtcenter_motionVstatic, mtperiphery_motionVstatic, mtall_motionVstatic
visualize(results.contrasts, sub, {'mtall_motionVstatic'})

%% save it (to draw ROIs, etc.)

hh = {'lh','rh'};
% get dimensions
hSize = get_surfsize(sub); % get hemi size
hSizeIdx=[1,hSize(1);hSize(1)+1,sum(hSize)];
% save mgz
origfile = fullfile(fullfile(projectDir, 'derivatives','freesurfer', sub, 'mri', 'orig.mgz'));
mgz = MRIread(origfile);

for hi=1:numel(hh)
    mgz.vol = results.contrasts.mtall_motionVstatic(hSizeIdx(hi,1):hSizeIdx(hi,2), 1);
    savePath = fullfile(projectDir, 'derivatives','freesurfer', sub, 'surf', sprintf('%s.motionloc.mgz',hh{hi}));
    MRIwrite(mgz, savePath);
end

disp('Manually draw the ROIs (labels and annot files required)')

%% after label / annot files are created:

surfaceROI = [];

for hi=1:numel(hh)
    % left hemi
    filename = sprintf('%s.mymotionloc.annot', hh{hi});
    [vertices, label, colortable] = read_annotation(fullfile(projectDir, 'derivatives','freesurfer', sub, 'label', filename));
    [idx,~] = find(strcmp(colortable.struct_names,'motionloc'));
    roiVal = colortable.table(idx,end);
    motion_idx = find(label==roiVal);
    motionmask = nan(length(vertices),1);
    motionmask(motion_idx) = 1;
    surfaceROI = [surfaceROI ; motionmask];
end

%%

v1ROI = [];

for hi=1:numel(hh)
    % left hemi
    filename = sprintf('%s.mymanualV1.annot', hh{hi});
    [vertices, label, colortable] = read_annotation(fullfile(projectDir, 'derivatives','freesurfer', sub, 'label', filename));
    [idx,~] = find(strcmp(colortable.struct_names,'manualV1'));
    roiVal = colortable.table(idx,end);
    v1_idx = find(label==roiVal);
    v1mask = nan(length(vertices),1);
    v1mask(v1_idx) = 1;
    v1ROI = [v1ROI ; v1mask];
end

% % just to check ROI in freesurfer (seems fine)
% for hi=1:numel(hh)
%     mgz.vol = surfaceROI(hSizeIdx(hi,1):hSizeIdx(hi,2), 1);
%     savePath = fullfile(projectDir, 'derivatives','freesurfer', sub, 'surf', sprintf('%s.motionlocSAVEDOUT.mgz',hh{hi}));
%     MRIwrite(mgz, savePath);
% end

%% Plot ROI I am using

vals = surfaceROI; %v1ROI; %
vals(isnan(vals)) = 0;
cvnlookup(sub,1,vals,[0 1],jet,0.5,[],[],{'overlayalpha',valsalpha,'roiname',{'Glasser2016'},'roicolor',{'b'},'drawroinames',1})

% setup color bar
title('Occipital Spherical Map of Manually Drawn MT+')
colormap(jet) % color map
caxis([0 1]) % range of the colorbar    
cb = colorbar();
cb.Ticks = 0:0.5:1;
ylabel(cb,'MT+ Mask','Rotation',270);
cb.Label.Position(1) = 5;
set(gca,'FontSize',20)

%% DRIFTING GRATING: load in datafiles for drifting grating

ses = 'ses-01';
stimdur_s = 3;
tr_s = 1;
run = 1:5;

% load_data (function in wpToolbox) - convert to .mgh if not already
datafiles = load_data(projectDir,projectName,dataStr,'.mgh',sub,ses,run);

designDir = fullfile(projectDir, 'designfiles');

% check if design files exist, if not convert the matrices from experiment
% into the appropriate format
matrices = format_desmats(projectDir, designDir, sub, ses, run, tr_s);


%% Visualize Drifting Grating Timeseries within MT

[nVertices, nTRs] = size((datafiles{1,1}));

for di=1:length(datafiles)
    baseline = mean(datafiles{di}, 2); % trying this
    timeseries = datafiles{di};
    timeseries_psc = (timeseries./baseline).*100;
    timeseries_psc = timeseries_psc-nanmean(timeseries_psc,'all');
    masked_fMRI = surfaceROI.* timeseries_psc; %datafiles{di,1};
    ave_fMRI = nanmean(masked_fMRI);
end

% you can see drift artifacts - must include tcompcor in GLM
figure
f = gcf;
f.Position = [-37 329 983 259];
plot(ave_fMRI, 'b','LineWidth', 2)
ylabel('% signal change')
xlabel('TR')
title('Mean signal in MT localizer')
ax = gca;
ax.FontSize = 14;

%% Trial triggered average

ttave(matrices,datafiles, surfaceROI, tr_s, stimdur_s)
%ttave(matrices,datafiles, v1ROI, tr_s, stimdur_s)

%% Run GLM for cardinal / oblique etc.

% per run
for rr=1:length(matrices)

    m={matrices{rr}};
    df={datafiles{rr}};

    % run glm (from GLMdenoise repo)
    results.(char(sprintf('allevents%s',num2str(rr)))) = GLMestimatemodel(m, df,...
        stimdur_s,tr_s,hrf_opt,[],resampling);
    
    % Script Part B) Description
    % Combines beta maps per condition and produces contrast maps:
    % average betas for cardinal motion - average betas for cardinal static
    % average betas for oblique motion - average betas for oblique static
    
    betamaps = results.(sprintf('allevents%s',num2str(rr))).modelmd{1,2};
    
    % average beta maps for cardinal moving - cardinal static
    contrasts.cardmVcards = mean(betamaps(:,1:4),2) - mean(betamaps(:,9:10),2);
    
    % average beta maps for oblique moving - oblique static
    contrasts.oblmVobls = mean(betamaps(:,5:8),2) - mean(betamaps(:,11:12),2);
    
    % average beta maps for all moving - all static
    contrasts.allmValls = mean(betamaps(:,1:8),2) - mean(betamaps(:,9:12),2);
    
    % average beta maps for all static - blank
    contrasts.allsVblank = mean(betamaps(:,9:12),2) - mean(betamaps(:,13),2);
    
    % average beta maps for all moving - blank
    contrasts.allmVblank = mean(betamaps(:,1:8),2) - mean(betamaps(:,13),2);
    
    contrasts.dir0mVdir0s = betamaps(:,1) - betamaps(:,9);
    contrasts.dir180mVdir180s = betamaps(:,3) - betamaps(:,9);
    contrasts.dir90mVdir90s = betamaps(:,2) - betamaps(:,10);
    contrasts.dir270mVdir270s = betamaps(:,4) - betamaps(:,10);
    contrasts.dir45mVdir45s = betamaps(:,5) - betamaps(:,11);
    contrasts.dir225mVdir225s = betamaps(:,7) - betamaps(:,11);
    contrasts.dir135mVdir135s = betamaps(:,6) - betamaps(:,12);
    contrasts.dir315mVdir315s = betamaps(:,8) - betamaps(:,12);

    results.(char(sprintf('contrasts%s',num2str(rr)))) = contrasts;
end

visualize(results.(char(sprintf('contrasts%s',num2str(rr)))), sub, {'allmValls'})

%% PART II: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RETINOTOPY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ses = 'ses-02';

% Description:
% This analyzes pRF results
% 1. load pRF data matrix and convert ecc to degrees

% load results
analyzePRFDir = fullfile(projectDir,'derivatives','analyzePRF', sub, ses, ...
    dataStr, modelType);
prfMainFile = strcat(sub, '_', ses, '_', modelType, '_results.mat');
prfresults = load(fullfile(analyzePRFDir,prfMainFile));
prfresults = prfresults.results;

% since results.ecc are in pixel values of stimulus - convert to stimulus radius
% range:
prfresults.ecc = ((prfresults.ecc).*stimradiusdeg)./100; % convert arbitrary units to degrees

% Take a look at the data dimensions

% prfresults.ang = polar angle values (0-360)
% prfresults.ecc = eccentricity values
% prfresults.rfsize = population receptive field size
% prfresults.R2 = variance explained (higher is better model fit)

% Check what size of one of these variables:
% Note that these variables are stored as fields within the struct variable
% "prfresults", so we call need to call the individual variable using "structName.fieldName" 
size(prfresults.ecc) % notice that "size" function returns you the size of each dimension of the input variable

% % % just to convert freesurfer surface to degrees (Angle)
for hi=1:numel(hh)
    tmp = prfresults.ang';
    if hi==1
        tmp(tmp>90 & tmp<270) = nan;
    end
    mgz.vol = tmp(hSizeIdx(hi,1):hSizeIdx(hi,2), 1);
    savePath = fullfile(projectDir, 'derivatives','freesurfer', sub, 'surf', sprintf('%s.angle_deg.mgz',hh{hi}));
    MRIwrite(mgz, savePath);
end

%% extract and visualize values of polar angle 

% disp('start the polar angle visualization')
% 
% % first convert everything to fsaverage
% metrics = {'ang', 'ecc', 'R2', 'rfsize'};
% setenv('SUBJECTS_DIR', [projectDir '/derivatives/freesurfer']); 
% hSize = get_surfsize(sub); % get hemi size
% 
% for mi=1:length(metrics)
%     vals = prfresults.(metrics{mi})'; % entire brain dimensions
%     hSizeIdx=[1,hSize(1);hSize(1)+1,numel(vals)];  
%     lhvals = double(vals(hSizeIdx(1,1):hSizeIdx(1,2)));  
%     rhvals = double(vals(hSizeIdx(2,1):hSizeIdx(2,2))); 
%     % convert to fsaverage space
%     fsaveragelh = cvntransfertosubject(sub,'fsaverage', lhvals, 'lh', 'nearest', 'orig', 'orig');
%     fsaveragerh = cvntransfertosubject(sub,'fsaverage', rhvals, 'rh', 'nearest', 'orig', 'orig');
%     fsaverageresults.(metrics{mi}) = [fsaveragelh ; fsaveragerh]';
% end


%% extract and visualize values of polar angle

vals = prfresults.R2'; %fsaverageresults.ang'; % choose variable to plot - needs to be a column
vals(prfresults.R2<.2) = nan; % select only certain valus of R2 (Nan for values lower than a threshold)
valsalpha = ~isnan(vals);
vals(isnan(vals)) = 0;
 
% show the results
cvnlookup(sub,1,vals,[0 100],jet,[],[],[],{'overlayalpha',valsalpha,'roiname',{'Glasser2016'},'roicolor',{'w'},'drawroinames',1})

% setup color bar
title('Occipital Spherical Map of R2 Values')
colormap(jet) % color map
caxis([0 100]) % range of the colorbar    
cb = colorbar();
cb.Ticks = 0:10:100;
ylabel(cb,'Variance Explained (%)','Rotation',270);
cb.Label.Position(1) = 5;
set(gca,'FontSize',20)


%% extract and visualize values of polar angle

vals = prfresults.ang'; %fsaverageresults.ang'; % choose variable to plot - needs to be a column
vals(prfresults.R2<.2) = nan; % select only certain valus of R2 (Nan for values lower than a threshold)
valsalpha = ~isnan(vals);
vals(isnan(vals)) = 0;
 
% show the results
cvnlookup(sub,1,vals,[0 360],hsv,[],[],[],{'overlayalpha',valsalpha,'roiname',{'Glasser2016'},'roicolor',{'b'},'drawroinames',1})

% setup color bar
title('Occipital Spherical Map of Polar Angle Values')
colormap(hsv) % color map
caxis([0 360]) % range of the colorbar    
cb = colorbar();
cb.Ticks = 0:45:360;
ylabel(cb,'Polar Angle (deg)','Rotation',270);
cb.Label.Position(1) = 5;
set(gca,'FontSize',20)



%% extract and visualize values of eccentricity

disp('start the eccentricity visualization')

vals = prfresults.ecc';
vals(prfresults.ecc>stimwidthdeg) = nan; % select eccentricity values < than the maximum eccentricity
vals(prfresults.ecc<100*eps) = nan; % this gets rid of the near zero values in medial surface
vals(prfresults.R2<.2) = nan; % select only certain valus of R2 (Nan for values lower than a threshold)
valsalpha = ~isnan(vals);
vals(isnan(vals)) = 0;
 
% show the results
cvnlookup(sub,1,vals,[0 stimradiusdeg],hsv,[],[],[],{'overlayalpha',valsalpha,'roiname',{'Glasser2016'},'roicolor',{'b'},'drawroinames',1})

% setup color bar
title('Occipital Spherical Map of Eccentricity Values')
colormap(flipud(jet)) % color map
caxis([0 stimradiusdeg]) % range of the colorbar    
cb = colorbar();
cb.Ticks = 0:5:stimradiusdeg;
ylabel(cb,'Eccentricity (deg of visual angle)','Rotation',270);
cb.Label.Position(1) = 5;
set(gca,'FontSize',20)

%%

% look inside MT/V1 to find polar angle & eccentricity representation
binsize = 45;
bintitles = {'0', '45', '90', '135', '180', '225', '270', '315'};
binlocations = [0 45 90 135 180 225 270 315];
binMasks = nan(length(binlocations),nVertices); % create masks per bin

ROI_angle = prfresults.ang.*surfaceROI'; %v1ROI'; %

for ri=1:numel(binlocations)
    binloc = binlocations(ri);
    if binloc == 0
        [~, vIndx] = find((ROI_angle <= binloc+binsize/2) | (ROI_angle >= 360-binsize/2)); 
    else
        [~, vIndx] = find((ROI_angle <= binloc+binsize/2) & (ROI_angle >= binloc-binsize/2)); 
    end
    length(vIndx)
    binMasks(ri, vIndx) = 1;
end

%%
figure
for ri=1:numel(binlocations)
    subplot(2,4,ri) 
    aveCon1 = nan(length(matrices),1); aveCon2 = nan(length(matrices),1);
    for rr=1:length(matrices)
        filterdata1 = results.(char(sprintf('contrasts%s',num2str(rr)))).cardmVcards'.*binMasks(ri, :);
        filterdata2 = results.(char(sprintf('contrasts%s',num2str(rr)))).oblmVobls'.*binMasks(ri, :);
        %filterdata1 = contrasts.cardmVcards'.*binMasks(ri, :);
        %filterdata2 = contrasts.oblmVobls'.*binMasks(ri, :);
        aveCon1(rr) = nanmean(filterdata1,2);
        aveCon2(rr) = nanmean(filterdata2,2);
    end
    allTrials_aveCon1 = mean(aveCon1);
    allTrials_aveCon2 = mean(aveCon2);
    allTrials_stdCon1 = std(aveCon1);
    allTrials_stdCon2 = std(aveCon2);
    y = [allTrials_aveCon1 allTrials_aveCon2]';
    h = bar(1, y);
    hold on
    xPos = [h(1).XData+h(1).XOffset; h(2).XData+h(2).XOffset];
    e = errorbar(xPos, y, [allTrials_stdCon1 allTrials_stdCon2]', '.k');
    title(sprintf('Wedge at %s deg', bintitles{ri}))
    ylim([-0.35 0.2]) %ylim([0 0.2])
    ax = gca;
    ax.XTick = [];
    %h.CData = [0 0 1; 1 0 0];
    set(h(1),'facecolor',[127, 191, 123]/255)
    set(h(2),'facecolor',[175, 141, 195]/255)
    if ri==1
        legend({'Cardinal Motion'; 'Oblique Motion'})
        ylabel('arbitrary units')
    elseif ri==5
        ylabel('arbitrary units')
    end
end
sgtitle('V1 Mean Motion Response per Polar Angle Locations')
f1 = gcf;
f1.Position = [12 402 1164 395];

%% TO DO: separate by each motion condition

figure
for ri=1:numel(binlocations)
    subplot(2,4,ri) 
    aveCon1 = nan(length(matrices),1); aveCon2 = nan(length(matrices),1);
    for rr=1:length(matrices)

        % location & direction match (radial outward)
        filterdata1_out = results.(char(sprintf('contrasts%s',num2str(rr)))).(char(sprintf('dir%smVdir%ss', bintitles{ri}, bintitles{ri})))'.*binMasks(ri, :);
        % location & direction opposite (radial2)
        inwarddir = binlocations(ri)-180;
        if inwarddir<0
            inwarddir = binlocations(ri)+180;
        end
        filterdata1_in = results.(char(sprintf('contrasts%s',num2str(rr)))).(char(sprintf('dir%smVdir%ss', num2str(inwarddir), num2str(inwarddir))))'.*binMasks(ri, :);
        filterdata1 = (filterdata1_out + filterdata1_in)./2;

        % tangential
        tangdir1 = median([binlocations(ri), inwarddir]);
        filterdata2_tang1 = results.(char(sprintf('contrasts%s',num2str(rr)))).(char(sprintf('dir%smVdir%ss', num2str(tangdir1), num2str(tangdir1))))'.*binMasks(ri, :);
        tangdir2 = tangdir1-180;
        if tangdir2<0
            tangdir2 = tangdir1+180;
        end
        filterdata2_tang2 = results.(char(sprintf('contrasts%s',num2str(rr)))).(char(sprintf('dir%smVdir%ss', num2str(tangdir2), num2str(tangdir2))))'.*binMasks(ri, :);
        filterdata2 = (filterdata2_tang1 + filterdata2_tang2)./2;

        aveCon1(rr) = nanmean(filterdata1,2);
        aveCon2(rr) = nanmean(filterdata2,2);
    end
    allTrials_aveCon1 = mean(aveCon1);
    allTrials_aveCon2 = mean(aveCon2);
    allTrials_stdCon1 = std(aveCon1);
    allTrials_stdCon2 = std(aveCon2);
    y = [allTrials_aveCon1 allTrials_aveCon2]';
    h = bar(1, y);
    hold on
    xPos = [h(1).XData+h(1).XOffset; h(2).XData+h(2).XOffset];
    e = errorbar(xPos, y, [allTrials_stdCon1 allTrials_stdCon2]', '.k');
    title(sprintf('Wedge at %s deg', bintitles{ri}))
    ylim([-.5 .4]) %ylim([0 0.25])
    ax = gca;
    ax.XTick = [];
    %h.CData = [0 0 1; 1 0 0];
    set(h(1),'facecolor',[146 197 222]/255)
    set(h(2),'facecolor',[202 0 32]/255)
    if ri==1
        legend({'Radial Motion'; 'Tangential Motion'})
        ylabel('arbitrary units')
    elseif ri==5
        ylabel('arbitrary units')
    end
end
sgtitle('V1 Mean Motion Response per Polar Angle Locations')
f1 = gcf;
f1.Position = [12 402 1164 395];




%% FIX LATER ~~~~~~~~~~~~

%% Visualize pRF/ecc by ROI

hh = {'lh','rh'}; % two hemispheres, left and right
nHemis = length(hh);

% check ROIfiles_Labeling.txt file in helper folder for the ROI index 
roiLabel = readtable(fullfile(projectDir, 'derivatives/freesurfer/fsaverage/atlasmgz/Glasser2016_ColorLUT.txt'));

% always set to fsaverage - just retrieved ROI indices (LH, RH)
% atlas can be changed inside function (currently Glasser 2016)
roiIdx = get_roi('fsaverage');

nRois = height(roiLabel(roiLabel.x_No_>=1,:)); 

% preallocate some cell structures to re-organize our data
rfSize = cell(nHemis,nRois); % pRF size
ecc = cell(nHemis,nRois);    % pRF eccentricity
R2 = cell(nHemis,nRois);     % R squared
polar = cell(nHemis,nRois);  % polar angle

% extract the number of verteces for left and right hemisphere (left -> row1; right -> row2)
%sub = 'fsaverage'; % take this out if in native space!!
%hSize = get_surfsize(sub); % get hemi size
%hSizeIdx=[1,hSize(1);hSize(1)+1,sum(hSize)];  % get hemi index where the maximum is the number of verteces with a functional activation (pRF)

for iH = 1:numel(hh) % for each hemisphere
    
    % get the temporary rfsize, ecc, and R2 from the orginal data (fsnative
    % space) fro one hemisphere
    tempRfsize = prfresults.rfsize(:,hSizeIdx(iH,1):hSizeIdx(iH,2));  
    tempEcc = prfresults.ecc(:,hSizeIdx(iH,1):hSizeIdx(iH,2));
    tempR2 = prfresults.R2(:,hSizeIdx(iH,1):hSizeIdx(iH,2));
    tempPolar = prfresults.ang(:,hSizeIdx(iH,1):hSizeIdx(iH,2));
    
    tempRfsize = tempRfsize';
    tempEcc = tempEcc';
    tempR2 = tempR2';
    tempPolar = tempPolar';
    % convert and update them into fsaverage space 
    tempRfsize = cvntransfertosubject(sub,'fsaverage', tempRfsize', hh{iH}, 'nearest', 'orig', 'orig');
    tempEcc = cvntransfertosubject(sub,'fsaverage', tempEcc', hh{iH}, 'nearest', 'orig', 'orig');
    tempR2 = cvntransfertosubject(sub,'fsaverage', tempR2', hh{iH}, 'nearest', 'orig', 'orig');
    tempPolar = cvntransfertosubject(sub,'fsaverage', tempPolar', hh{iH}, 'nearest', 'orig', 'orig');    
        
    for iROI = 1:nRois
        
        % create an index with the verteces that belong to a specific ROI
        indx = roiIdx(:,iH)==iROI;
        
        rfSize{iH,iROI} = tempRfsize(indx);
        ecc{iH,iROI} = tempEcc(indx);
        R2{iH,iROI} = tempR2(indx);
        polar{iH,iROI} = tempPolar(indx);
        
    end
end

%%

% used for thresholding
R2Min = 2; % set a minimum value for R squared
eccMax = stimradiusdeg; % set a maximum eccentricity

% plot relevant ROIs:
chosenROIs = {'V1','V2','V3','MT','LO1','LO2'} ; %'V4','MT','MST'};
figure

cmap = jet();
color_incr = floor(length(cmap)/length(chosenROIs)*0.8);
color_init = 0;

% all data is concatentated across hemis: [LH,RH]
for iR=1:length(chosenROIs)
    color = cmap(color_incr*iR,:); % colorvalues{iR};
    currROI = chosenROIs{iR};
    roiIdx = roiLabel(strcmp(roiLabel.LabelName_,currROI),:).x_No_;

    ecc_unfiltered = [ecc{1,roiIdx}; ecc{2,roiIdx}];
    rf_unfiltered = [rfSize{1,roiIdx}; rfSize{2,roiIdx}];
    R2_unfiltered = [R2{1,roiIdx}; R2{2,roiIdx}];
    polar_unfiltered = [polar{1,roiIdx}; polar{2,roiIdx}];

    % filter by eccentricity thresh
    [I1, ~] = find(ecc_unfiltered>eccMax*2);
    ecc_filtered = ecc_unfiltered; ecc_filtered(I1,:)=NaN;
    rf_filtered = rf_unfiltered; rf_filtered(I1,:)=NaN;
    R2_filtered = R2_unfiltered; R2_filtered(I1,:)=NaN;
    polar_filtered = polar_unfiltered; polar_filtered(I1,:)=NaN;
    
    % filter by R2 thresh
    [I2, ~] = find(R2_unfiltered<R2Min);
    ecc_filtered(I2,:)=NaN;
    rf_filtered(I2,:)=NaN;
    R2_filtered(I2,:)=NaN;
    polar_filtered(I2,:)=NaN;

    % fitting a line, normally "lsline" does the same job but this makes
    % the figure easier to manipulate
    idx = ~isnan(ecc_filtered);
    s = scatter(ecc_filtered(idx),rf_filtered(idx), 25, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.5);
    s.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on

    %find slope
    p = fitlm(ecc_filtered(idx),rf_filtered(idx),'Intercept',false);
    fitest = p.Coefficients.Estimate;
    x = linspace(0,eccMax); % make some x value
    y = fitest.*x;
    plot(x,y,'-','LineWidth',4,'Color',color); % plot the fitted line for this ROI

end

xlim([0 eccMax])
legend(chosenROIs)
xlabel('eccentricity (degrees)')
ylabel('pRF size (degrees)')
title('Relation between pRF Size & Eccentricity Across ROIs')
ax = gca;
ax.FontSize = 14;



%% Helper Function: COPIED FROM bidsAnalyzePRF_update.m

function [averageScans, stimwidth, opt] = getPRFOpts(prfOptsPath)

if ~exist('prfOptsPath', 'var') || isempty(prfOptsPath)
    prfOptsPath = prfOptsMakeDefaultFile;
end

json = jsondecode(fileread(prfOptsPath));


if isfield(json, 'averageScans'), averageScans = json.averageScans;
else, averageScans = []; end

if isfield(json, 'stimwidth'), stimwidth = json.stimwidth;
else, error('Stim width not specified in options file'); end

if isfield(json, 'opt'), opt = json.opt; else, opt = []; end

end

function pth = prfOptsMakeDefaultFile()
    % see analyzePRF for descriptions of optional input

    % average scans with identical stimuli
    json.averageScans = [];  %
    json.stimwidth    = 24.2;  % degrees

    % other opts
    json.opt.vxs            = [];
    json.opt.wantglmdenoise = [];
    json.opt.hrf            = [];
    json.opt.maxpolydeg     = [];
    json.opt.numperjob      = [];
    json.opt.xvalmode       = [];
    json.opt.seedmode       = [];
    json.opt.maxiter        = [];
    json.opt.display        = 'off';
    json.opt.typicalgain    = [];


    pth = fullfile(tempdir, 'prfOpts.json');
    savejson('', json, 'FileName', pth);
end