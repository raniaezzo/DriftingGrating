function [task, frameCounter, vbl] = my_stim(my_key, scr,const,expDes, task, frameCounter, trialID, vbl)

    trialType = expDes.trialMat(trialID,2);
    
    if trialType~=0   
        %tic
        if trialType==1  
            angle = const.maporientation(expDes.trialMat(trialID,3)); % 3rd column is orientation
            sprintf('ORIENTATION %i', expDes.trialMat(trialID,3))
        elseif trialType==2
            angle = const.mapdirection(expDes.trialMat(trialID,4)); % 4th column is direction
            sprintf('MOVING %i', expDes.trialMat(trialID,4))
        end

        waitframes = 1;
        waitduration = waitframes * scr.ifi;
        shiftperframe= const.stimSpeed_cps * const.stimSpeed_ppc * waitduration;
     

        movieDurationSecs=expDes.stimDur_s;   % Abort after 3 seconds.
        i = const.phaseLine(trialID);
        
        xoffset = mod(i*shiftperframe,const.stimSpeed_ppc);
        visiblesize = const.visiblesize;
        srcRect=[xoffset 0 xoffset + visiblesize visiblesize];
        gratingtex = const.gratingtex;

        xDist = 0; yDist = 0;
        xDist = scr.windCenter_px(1)+xDist-(visiblesize/2); % center + (+- distance added in pixels)
        yDist = scr.windCenter_px(2)+yDist-(visiblesize/2);  % check with -(vis part.. 
        dstRect=[xDist yDist visiblesize+xDist visiblesize+yDist];


        
        %vbl=Screen('Flip', const.window);
        vblendtime = vbl + movieDurationSecs;
        %i=0;
        
        if const.EL_mode, Eyelink('message', 'STIMULUS ONSET'); end
        % Animationloop:
        while (vbl < vblendtime)
        
        
            % Draw grating texture, rotated by "angle":
            Screen('DrawTexture', const.window, gratingtex, srcRect, dstRect, angle, ...
                [], [], [], [], []); %, propertiesMat');


            % Draw stimuli here, better at the start of the drawing loop
            fixColor = const.black;
            my_fixation(scr,const,fixColor)

            %vbl = Screen('Flip', const.window, vbl + (waitframes - 0.5) * scr.ifi);

            % pasted below ~~~~~~~~~
            %while (trialStart <= vbl) && (vbl <= trialStart + 4.5) && (runFinish~=1) %4.5 sec

            Screen('DrawingFinished',const.window); % small ptb optimisation

            vbl = Screen('Flip',const.window, vbl + (waitframes - 0.5) * scr.ifi);
        end
        
    end


end
