%% define path
clearvars; clc; close all;
subject = 'sub-0248';
bidsDir = '/Volumes/Vision/MRI/recon-bank';
%bidsDir = '~/Documents/MRI/bigbids';

%% tmp
view_fv(subject,bidsDir,'mt+2')
view_fv(subject,bidsDir,'mt+2','cd2','oppo3','prfvista_mov/eccen','prfvista_mov/angle_adj','myelin/SmoothedMyelinMap_BC');

%% compare 2D vs 3D 
view_fv(subject,bidsDir,'mt+2','cd2');
% view_fv(sub,bidsDir,'mt+12','cdavg');

%% compare four different motion condition 
view_fv(subject,bidsDir,'out','in','cw','ccw');

%% compare upper and lower meridian
view_fv(subject,bidsDir,'upper','lower','mt+2');

%% compare new and old mt+ localizer 
view_fv(subject,bidsDir,'mt+1','mt+2');

%% compare new and old cd 
view_fv(subject,bidsDir,'cd1','cd2');

%% compare new and old transparent motion 
view_fv(subjectDir,resultsDir,'oppo2','oppo3');

%% compare transparent motion vs. 3D motion 
view_fv(subject,bidsDir,'oppo3','cd2');
% view_fv(sub,bidsDir,'oppo3','cdavg');

%% myelin map 
view_fv(subject,bidsDir,'myelin/SmoothedMyelinMap_BC','motion_base/mt+2');
%view_fv(sub,bidsDir,'SmoothedMyelinMap','mt+2');
%view_fv(sub,bidsDir,'MyelinMap_BC','mt+2');
%view_fv(sub,bidsDir,'MyelinMap','mt+2');

%% prf
%view_fv(subject,bidsDir,'angle_adj','eccen','sigma');
cmd = view_fv(subject,bidsDir,'prfvista_mov/eccen');
%% compare how many runs to use
view_fv(subject,bidsDir,'run4','run3','run2','run1'); % which runs
view_fv(subject,bidsDir,'4run','3run','2run','run1'); % how many runs

%% hand vs. mt+
view_fv(subject,bidsDir,'hand','mt+2');
