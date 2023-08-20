function [task, frameCounter, vbl] = my_blank(my_key, scr,const, task, frameCounter, interval, vbl)

try
    waitframes = 1;
    %vbl = Screen('Flip',const.window);
    vblendtime = vbl + interval;

    % Blank period
    while vbl <= vblendtime  

        if task(frameCounter,1)==1
            fixColor = const.lightgray;
        else
            fixColor = const.black;
        end

        % draw stimuli here, better at the start of the drawing loop
        Screen('DrawTexture', const.window, const.pinknoiseTex); % just added
        my_fixation(scr,const,fixColor)
        Screen('DrawingFinished',const.window); % small ptb optimisation
        vbl = Screen('Flip',const.window, vbl + (waitframes - 0.5) * scr.ifi);

        % check for keyboard input
        [keyIsDown, ~, keyCode] = KbCheck;
        if ~keyIsDown
            [keyIsDown, ~, keyCode] = KbCheck(my_key.suppResponseID);
        end
        if keyIsDown && keyCode(my_key.escape)
            ShowCursor; sca; return
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
    
catch
    vpixxShutdown(const); 
end

end