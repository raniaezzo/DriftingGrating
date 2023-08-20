function [task, frameCounter, vbl] = my_stim(my_key, scr,const,expDes, task, frameCounter, trialID, vbl)

try
    trialType = expDes.trialMat(trialID,2); % 0=baseline, 1=static, 2=motion

    if trialType~=0   

        %tic
        if trialType==1 
            angle = const.maporientation(expDes.trialMat(trialID,3)); % 3rd column is orientation
        elseif trialType==2
            angle = const.mapdirection(expDes.trialMat(trialID,4)); % 4th column is direction
        end

        movieDurationSecs=expDes.stimDur_s;   % Abort after 3 seconds.
        i = const.phaseLine(trialID);

        % included: ~~~~~

        grating_halfw = const.grating_halfw;
        visiblesize = const.visiblesize;
        gratingtex = const.gratingtex;
        maskOutertex = const.maskOutertex;
        maskInnertex = const.maskInnertex;

        % eccentricity = 0
        xDist = 0; yDist = 0;
        xDist = scr.windCenter_px(1)+xDist-(visiblesize/2); % center + (+- distance added in pixels)
        yDist = scr.windCenter_px(2)+yDist-(visiblesize/2);  % check with -(vis part.. 
        dstRect=[xDist yDist visiblesize+xDist visiblesize+yDist];

        % % properties matrix
        % propertiesMat = [const.phaseLine(1), const.stimSF_cpp, ...
        %     const.gaussianSigma, const.contrast, const.aspectRatio, 0, 0, 0];

        waitframes = 1;
        waitduration = waitframes * scr.ifi;
        shiftperframe= const.stimSpeed_cps * const.stimSpeed_ppc * waitduration;

        %vbl=Screen('Flip', const.window);
        vblendtime = vbl + movieDurationSecs;
        %i=0;

        if const.EL_mode, Eyelink('message', 'STIMULUS ONSET'); end
        % Animationloop:
        while (vbl < vblendtime)

            if task(frameCounter,1)==1
                fixColor = const.lightgray;
            else
                fixColor = const.black;
            end

            xoffset = mod(i*shiftperframe,const.stimSpeed_ppc);

            if trialType==2   
                i=i+1;
            end

            srcRect=[xoffset 0 xoffset + visiblesize visiblesize];

            % Set the right blend function for drawing the gabors
            %Screen('BlendFunction', const.window, 'GL_ONE', 'GL_ZERO');

            % Draw grating texture, rotated by "angle":
            Screen('DrawTexture', const.window, gratingtex, srcRect, dstRect, angle, ...
                [], [], [], [], []); %, propertiesMat');

            % outer and inner masks
            Screen('DrawTexture', const.window, maskOutertex, [], [], []); %[0 0 scr.windX_px scr.windY_px], []);
            Screen('DrawTexture', const.window, maskInnertex, [], dstRect, []);

            % Change the blend function to draw an antialiased fixation point
            % in the centre of the array
            %Screen('BlendFunction', const.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

            % Draw stimuli here, better at the start of the drawing loop
            my_fixation(scr,const,fixColor)

            %vbl = Screen('Flip', const.window, vbl + (waitframes - 0.5) * scr.ifi);

            % pasted below ~~~~~~~~~
            %while (trialStart <= vbl) && (vbl <= trialStart + 4.5) && (runFinish~=1) %4.5 sec

            Screen('DrawingFinished',const.window); % small ptb optimisation

            vbl = Screen('Flip',const.window, vbl + (waitframes - 0.5) * scr.ifi);

            % check for keyboard input
            [keyIsDown, ~, keyCode] = KbCheck(my_key.keyboardID);
            if ~keyIsDown
                [keyIsDown, ~, keyCode] = KbCheck(my_key.suppResponseID);
            end
            if keyIsDown && keyCode(my_key.escape)
                ShowCursor; sca; clear mex; clear fun; return
            elseif keyIsDown && ~keyCode(my_key.escape) && ~(keyCode(my_key.Trigger) || keyCode(34))
                task(frameCounter,2) = 1;   
            end

            % vpixx
            if const.vpixx == 1
                Datapixx('RegWrRd');
                kbcheck = dec2bin(Datapixx('GetDinValues'));
                if kbcheck(end) == '1' || kbcheck(end-3) == '1' || kbcheck(end-2) == '1' || kbcheck(end-1) == '1'|| kbcheck(end-4) == '1'|| kbcheck(end-5) == '1'|| kbcheck(end-6) == '1'|| kbcheck(end-7) == '1'|| kbcheck(end-8) == '1'|| kbcheck(end-9) == '1'
                %if kbcheck(end) == '1' || kbcheck(end-3) == '1' || kbcheck(end-2) == '1' || kbcheck(end-1) == '1'
                    task(frameCounter,2) = 1;  
                end
            end

            FlushEvents('KeyDown');
            frameCounter=frameCounter+1;

        end

        if const.EL_mode, Eyelink('message', 'STIMULUS OFFSET'); end
        rep = 1; % only run iti (not a longer blank interval)
        %toc
    else
        %rep = 3;
        rep = (expDes.stimDur_s+expDes.itiDur_s)/expDes.itiDur_s;
    end

    % shared code for baseline & also intertrial interval
    if const.EL_mode, Eyelink('message', 'BLANK (2 s) START'); end
    [task, frameCounter, vbl] = my_blank(my_key, scr,const, task, frameCounter, expDes.itiDur_s*rep, vbl);
    if const.EL_mode, Eyelink('message', 'BLANK (2 s) END'); end  
    
catch
    vpixxShutdown(const); return
end

end
