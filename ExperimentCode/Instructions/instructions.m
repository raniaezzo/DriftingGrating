function keyCode = instructions(scr,const,my_key,text)
% ----------------------------------------------------------------------
% instructions(scr,const,my_key,text,button)
% ----------------------------------------------------------------------
% Goal of the function :
% Display instructions write in a specified matrix.
% ----------------------------------------------------------------------
% Input(s) :
% scr : main window pointer.
% const : struct containing all the constant configurations.
% text : library of the type {}.
% ----------------------------------------------------------------------
% Output(s):
% (none)
% ----------------------------------------------------------------------

x_mid = scr.windCenter_px(1);
y_mid = scr.windCenter_px(2);

%while KbCheck(-1); end
while KbCheck(my_key.keyboardID); end
KbName('UnifyKeyNames');

push_button = 0;
while ~push_button
    
    Screen('Preference', 'TextAntiAliasing',1);
    Screen('TextSize',const.window, const.text_size);
    Screen ('TextFont', const.window, const.text_font);
    Screen('FillRect', const.window, const.gray);
    
    sizeT = size(text);
    lines = sizeT(1)+2;
    bound = Screen('TextBounds',const.window,text{1,:});
    espace = ((const.text_size)*1.50);
    first_line = y_mid - ((round(lines/2))*espace);
    
    addi = 0;
    for t_lines = 1:sizeT(1)
        Screen('DrawTexture', const.window, const.pinknoiseTex); % just added
        Screen('DrawText',const.window,text{t_lines,:},x_mid-bound(3)/2,first_line+addi*espace, const.white);
        addi = addi+1;
    end
    addi = addi+2;
    my_fixation(scr,const,const.black)
    Screen('Flip',const.window);
    
    % wait for trigger with keyboard (keep with VPIXX for ability to quit)
    [ keyIsDown, ~, keyCode ] = KbCheck(my_key.keyboardID); % KbCheck(-1);
    if keyIsDown
        disp('KEYCODE IS')
        disp(keyCode)
        if keyCode(my_key.Trigger) || keyCode(34)
            disp('TRIGGERED!')
            push_button=1;
        elseif keyCode(my_key.escape) && ~const.expStart
            return
        end
    end
    
    if const.vpixx == 1
        % wait for trigger from datapixx
        Datapixx('RegWrRd');
        init_check = dec2bin(Datapixx('GetDinValues'));
        trigger_state = init_check(14);
        while 1
            Datapixx('RegWrRd');
            regcheck = dec2bin(Datapixx('GetDinValues'));
            if regcheck(14) ~= trigger_state
                push_button=1;
                disp('TRIGGERED!')
                break;
            end
        end
    end

end

% keyPressed = find(keyCode);

end
