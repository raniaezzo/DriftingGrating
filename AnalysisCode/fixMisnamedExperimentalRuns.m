% fix incorrect subject and run number
% assume files have been renamed



%% to rename subject

% path = '/Volumes/Vision/UsersShare/Rania/Project_dg/experimentalOutput/da/02';
% runname = 'Run1';
% subjectid = '02'; misNamedid = '00';
% load(fullfile(path, runname, sprintf('S%s_const_file_%s.mat', subjectid, runname)));

% const.subjID = strrep(const.subjID, misNamedid, subjectid);
% const.subjDir = strrep(const.subjDir, misNamedid, subjectid);
% const.runLog  = strrep(const.runLog, misNamedid, subjectid);
% const.runDir  = strrep(const.runDir, misNamedid, subjectid);
% const.scr_fileDat  = strrep(const.scr_fileDat, misNamedid, subjectid);
% const.scr_fileMat  = strrep(const.scr_fileMat, misNamedid, subjectid);
% const.const_fileDat  = strrep(const.const_fileDat, misNamedid, subjectid);
% const.const_fileMat  = strrep(const.const_fileMat, misNamedid, subjectid);
% const.expRes_fileCsv  = strrep(const.expRes_fileCsv, misNamedid, subjectid);
% const.design_fileMat  = strrep(const.design_fileMat, misNamedid, subjectid);
% const.responses_fileMat  = strrep(const.responses_fileMat, misNamedid, subjectid);
% const.eyeDataDir = strrep(const.eyeDataDir, misNamedid, subjectid);
% save(fullfile(path, runname, sprintf('S%s_const_file_%s.mat', subjectid, runname)), 'const');

%% to rename run 
path = '/Volumes/Vision/UsersShare/Rania/Project_dg/experimentalOutput/da/02';

runname = 'Run8'; runNumer = strsplit(runname, 'Run'); runNumer = runNumer{2};
misNamedrun = 'Run7';

load(fullfile(path, runname, sprintf('S%s_const_file_%s.mat', subjectid, runname)));

const.run = str2double(runNumer);
const.runDir  = strrep(const.runDir, misNamedrun, runname);
const.scr_fileDat  = strrep(const.scr_fileDat, misNamedrun, runname);
const.scr_fileMat  = strrep(const.scr_fileMat, misNamedrun, runname);
const.const_fileDat  = strrep(const.const_fileDat, misNamedrun, runname);
const.const_fileMat  = strrep(const.const_fileMat, misNamedrun, runname);
const.expRes_fileCsv  = strrep(const.expRes_fileCsv, misNamedrun, runname);
const.design_fileMat  = strrep(const.design_fileMat, misNamedrun, runname);
const.responses_fileMat  = strrep(const.responses_fileMat, misNamedrun, runname);
const.eyeDataDir = strrep(const.eyeDataDir, misNamedrun, runname);

save(fullfile(path, runname, sprintf('S%s_const_file_%s.mat', subjectid, runname)), 'const');


