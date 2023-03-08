function attentiontask = my_task(expDes, scr)

    % generate random fixation change for entire run (gamma)

    % define parameters in terms of frames, not secs
    changeduration_frames = round(0.5/scr.ifi);
    shape = 0.1 / scr.ifi;
    scale = 0.4 / scr.ifi; % smaller = more frequent events
    % 0.2 and 0.5 are good for ~8x per sec
    
    events = zeros(expDes.totalframes,1);
    eventend = 0; eventstart=0;

    while 1
        % define start of event
        eventstart = eventend + round(random('gamma',shape,scale));

        % define end of event
        eventend = eventstart+changeduration_frames;

        % add ones for event occurance
        if eventend <= length(events)
            events(eventstart:eventend,1) = 1;
        else
            break
        end

    end

    responses = zeros(expDes.totalframes,1);
    attentiontask = [events, responses];
    
end


