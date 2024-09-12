function [task, frameCounter, vbl] = my_stim(my_key, scr,const,expDes, task, frameCounter, trialID, vbl)
phase=0;
try
    trialType = expDes.trialMat(trialID,2); % 0=baseline, 1=static, 2=motion

    if trialType~=0   

        %tic
        if trialType==1 
            angle = const.maporientation(expDes.trialMat(trialID,3)); % 3rd column is orientation
            sprintf('ORIENTATION %i', expDes.trialMat(trialID,3))
        elseif trialType==2
            angle = const.mapdirection(expDes.trialMat(trialID,4)); % 4th column is direction
            sprintf('MOVING %i', expDes.trialMat(trialID,4))
        end

        movieDurationSecs=expDes.stimDur_s;   % Abort after 3 seconds.
        i = const.phaseLine(trialID);

        % included: ~~~~~

        grating_halfw = const.grating_halfw;
        visiblesize = const.visiblesize;
        
        gratingtex = const.gratingtex;
        
        if strcmp(const.expType, 'dgl')
            if trialType==1 % static
                phaseJump = 0; 
                orientation = const.dglmaporientation(expDes.trialMat(trialID,3));
            elseif trialType==2 % moving
                phaseJump = 25; % always this way (unlike da)
                orientation = const.dglmapdirection(expDes.trialMat(trialID,4));
            end
        elseif strcmp(const.expType, 'da')
            if trialType==1 % static
                phaseJump = 0; 
            elseif trialType==2 % moving
                if ismember(expDes.trialMat(trialID,4), [0, 45, 270, 315]) % if moving clockwise (includes inward)
                    phaseSign = 1;
                elseif ismember(expDes.trialMat(trialID,4), [90, 180, 135, 225]) % if movibg counterclockwise (includes outward)
                    phaseSign = -1;
                end
                phaseJump = phaseSign*25; % this is only used for 'da'
            end
            virtualSize = const.grating_halfw*2;
            frequency = const.stimSF_cpp; %0.02; % cycles/pixel (0.06)
            middleRadius = virtualSize/2;
            middlePerimeter = 2*pi*middleRadius; % pixels
            % has to scale with eccentricity (half way point)
            radialFrequency = (frequency*middlePerimeter / (2*pi))*.5;
            % radial orientation (pinwheel)
            if (trialType==1 && (expDes.trialMat(trialID,3) == 90)) || ...
               (trialType==2 && (expDes.trialMat(trialID,4) == 0)) || ...
               (trialType==2 && (expDes.trialMat(trialID,4) == 180))
                circularFrequency = 0;
            % tangential orientation (annulus)
            elseif (trialType==1 && (expDes.trialMat(trialID,3) == 0)) || ...
                (trialType==2 && (expDes.trialMat(trialID,4) == 90)) || ...
                (trialType==2 && (expDes.trialMat(trialID,4) == 270))
                radialFrequency = 0;
                circularFrequency = frequency;
            % oblique (spiral)
            elseif (trialType==1 && (expDes.trialMat(trialID,3) == 135)) || ...
               (trialType==2 && (expDes.trialMat(trialID,4) == 45)) || ...
               (trialType==2 && (expDes.trialMat(trialID,4) == 225))
                circularFrequency = -frequency; 
            elseif (trialType==1 && (expDes.trialMat(trialID,3) == 45)) || ...
               (trialType==2 && (expDes.trialMat(trialID,4) == 135)) || ...
               (trialType==2 && (expDes.trialMat(trialID,4) == 315))
                circularFrequency = frequency; 
            end
        end
        
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
            
            % just added
            if strcmp(const.expType, 'dgl')
                Screen('DrawTexture', const.window, const.pinknoiseTex, [], [], []);
            end

            if task(frameCounter,1)==1
                fixColor = const.lightgray;
            else
                fixColor = const.black;
            end

            xoffset = mod(i*shiftperframe,const.stimSpeed_ppc);

            if trialType==2 %&& strcmp(const.expType, 'dg')
                i=i+1;
            end

            srcRect=[xoffset 0 xoffset + visiblesize visiblesize];

            % Set the right blend function for drawing the gabors
            %Screen('BlendFunction', const.window, 'GL_ONE', 'GL_ZERO');

            if strcmp(const.expType, 'dg')
                % Draw grating texture, rotated by "angle":
                Screen('DrawTexture', const.window, gratingtex, srcRect, dstRect, angle, ...
                    [], [], [], [], []); %, propertiesMat');
            elseif strcmp(const.expType, 'da')
                phase = phase + phaseJump;
                angle = 0;
                Screen('DrawTexture', const.window, gratingtex, [], dstRect,...
                    angle, [], [], [0.5 0.5 0.5 1], [], [],...
                    [phase, radialFrequency, 0.5, -1, circularFrequency, 0, 0, 0]);
                    % for LAST [] AUXILIARY PARAMS:
                    % input3: contrast
                    % input4: < 0 is a sinusoid.
                    % input5: circular frequency
            elseif strcmp(const.expType, 'dgl')
                phase = phase + phaseJump; % phase = 0; 
                baseColor = [0.5 0.5 0.5 1];
                frequency = 0.02; contrast = 0.5;
                
                Screen('BlendFunction', const.window, 'GL_ONE', 'GL_ZERO');
                % UVM position (angle = 0); 
                Screen('DrawTexture', const.window, gratingtex, [], dstRect, ...
                    0, [], [], baseColor, [], [], ...
                    [phase, frequency, contrast, orientation]);
                
                % UVM position (angle = 0); 
                Screen('DrawTexture', const.window, gratingtex, [], dstRect, ...
                    45, [], [], baseColor, [], [], ...
                    [phase, frequency, contrast, orientation]);
                
                % UVM position (angle = 0); 
                Screen('DrawTexture', const.window, gratingtex, [], dstRect, ...
                    90, [], [], baseColor, [], [], ...
                    [phase, frequency, contrast, orientation]);
                
                % UVM position (angle = 0); 
                Screen('DrawTexture', const.window, gratingtex, [], dstRect, ...
                    135, [], [], baseColor, [], [], ...
                    [phase, frequency, contrast, orientation]);
                
                % UVM position (angle = 0); 
                Screen('DrawTexture', const.window, gratingtex, [], dstRect, ...
                    180, [], [], baseColor, [], [], ...
                    [phase, frequency, contrast, orientation]);
                
                % UVM position (angle = 0); 
                Screen('DrawTexture', const.window, gratingtex, [], dstRect, ...
                    225, [], [], baseColor, [], [], ...
                    [phase, frequency, contrast, orientation]);
                
                % UVM position (angle = 0); 
                Screen('DrawTexture', const.window, gratingtex, [], dstRect, ...
                    270, [], [], baseColor, [], [], ...
                    [phase, frequency, contrast, orientation]);
                
                % UVM position (angle = 0); 
                Screen('DrawTexture', const.window, gratingtex, [], dstRect, ...
                    315, [], [], baseColor, [], [], ...
                    [phase, frequency, contrast, orientation]);
            end
            
            Screen('BlendFunction', const.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            % outer and inner masks
            Screen('DrawTexture', const.window, maskOutertex, [], [], []); %[0 0 scr.windX_px scr.windY_px], []);
            Screen('BlendFunction', const.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
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
