function overDone(const, responses, trial_onsets)
% ----------------------------------------------------------------------
% overDone
% ----------------------------------------------------------------------
% Goal of the function :
% Close screen, listen keyboard and save duration of the experiment
% ----------------------------------------------------------------------
% Input(s) :
% none
% ----------------------------------------------------------------------
% Output(s):
% none
% ----------------------------------------------------------------------

fid = fopen(const.runLog, 'a');
fprintf(fid, '\n')
fprintf(fid, 'Run%i', const.run)

% .mat file
save(const.responses_fileMat,'responses');

timestamp_path = strrep(const.responses_fileMat,'responses', 'timestamps');
save(timestamp_path,'trial_onsets');

%figure
%plot(responses(1:length(responses)/4))

ListenChar(1);
%WaitSecs(2.0);

ShowCursor;
Screen('CloseAll');
sca;
timeDur=toc/60;
fprintf(1,'\nTotal time : %2.0f min.\n\n',timeDur);

%PsychPortAudio('Stop', const.pahandle);
%PsychPortAudio('Close', const.pahandle); %added

clear mex;
clear fun;

end

