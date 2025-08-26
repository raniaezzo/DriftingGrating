function output = detectSacc(d,vel,samplerate)

    % adapted from FUNCTION microsacc.m
    %  Detection of monocular candidates for microsaccades;
    %  Please cite: Engbert, R., & Mergenthaler, K. (2006) Microsaccades 
    %  are triggered by low retinal image slip. Proceedings of the National 
    %  Academy of Sciences of the United States of America, 103: 7192-7197.
    
    % and microsaccMerge.m from Ezzo, R. Song B., Rokers, B., & Carrasco M.
    % (2025)
    VFAC=[6,6]; % velocity threshold for X, Y to constitute as a microsaccade
    MINDUR=6; % minimum duration (ms) to constitute as a microsaccade 
    mergeInterval=15; % minimum duration (ms) to constitute as separate microsaccades; otherwise merged as one 
    threshold_AM=[.05 1]; 
    
    % compute threshold
    msdx = sqrt( median(vel(:,1).^2) - (median(vel(:,1)))^2 );
    msdy = sqrt( median(vel(:,2).^2) - (median(vel(:,2)))^2 );

    radiusx = VFAC(1)*msdx;
    radiusy = VFAC(2)*msdy;
    radius = [radiusx radiusy];
    
    % compute test criterion: ellipse equation
    test = (vel(:,1)./radiusx).^2 + (vel(:,2)./radiusy).^2;
    indx = find(test>1);
    
    % determine saccades
    N = length(indx); 
    sac = [];
    nsac = 0;
    dur = 1;
    a = 1;
    k = 1;
    while k<N
        if indx(k+1)-indx(k)==1
            dur = dur + 1;
        else
            if dur>=MINDUR*samplerate/1000
                nsac = nsac + 1;
                b = k;
                sac(nsac,:) = [indx(a) indx(b)];
            end
            a = k+1;
            dur = 1;
        end
        k = k + 1;
    end

    % merge saccades
    if ~isempty(sac)
        msac = sac(1,:);    % merged saccade matrix
        s    = 1;           % index of saccades in sac
        sss  = 1;           % boolean for still same saccade
        nsac = 1;           % number of saccades after merge
        while s<size(sac,1)
            if ~sss
                nsac = nsac + 1;
                msac(nsac,:) = sac(s,:);
            end
            if sac(s+1,1)-sac(s,2) <= mergeInterval*samplerate/1000
                msac(nsac,2) = sac(s+1,2);
                sss = 1;
            else
                sss = 0;
            end
            s = s+1;
        end
        if ~sss
            nsac = nsac + 1;
            msac(nsac,:) = sac(s,:);
        end
    else
        msac = [];
        nsac = 0;
    end

    % compute peak velocity, horizonal and vertical components
    for s=1:nsac
        % onset and offset
        a = msac(s,1); 
        b = msac(s,2); 
        % saccade peak velocity (vpeak)
        vpeak = max( sqrt( vel(a:b,1).^2 + vel(a:b,2).^2 ) );
        msac(s,3) = vpeak;
        % saccade vector (dx,dy)
        dx = d(b,1)-d(a,1); 
        dy = d(b,2)-d(a,2); 
        msac(s,4) = dx;
        msac(s,5) = dy;
        % saccade amplitude (dX,dY)
        i = msac(s,1):msac(s,2);
        [minx, ix1] = min(d(i,1));
        [maxx, ix2] = max(d(i,1));
        [miny, iy1] = min(d(i,2));
        [maxy, iy2] = max(d(i,2));
        dX = sign(ix2-ix1)*(maxx-minx);
        dY = sign(iy2-iy1)*(maxy-miny);
        amplitude=sqrt((dX.^2)+(dY.^2)); % compute amplitude from X, Y components
        msac(s,6:7) = [dX dY];
        msac(s,8) = amplitude;
        
    end

    % filtered_msac includes the filter - in case I want to do msacc analysis
        % later ('msac' has no filter)
    if ~isempty(msac)
        order= msac(:,8) > threshold_AM(1) & msac(:,8) < threshold_AM(2);
        filtered_msac= msac(order,:);
    else
        msac=[]; filtered_msac = [];
    end

    output = {msac, filtered_msac};

end