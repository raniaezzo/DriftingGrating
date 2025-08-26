function metaDataStruct = extractmsgMetadata(filename, subjectname)

    % find lines that contain "RATE" in msg file
    A = regexp(fileread(filename),'\n','split');
    
    whichline = find(contains(A,'RATE'));
    RATELINE = split(char(A(whichline(1))));
    samplingRate = str2double(RATELINE{5});
    
    
    whichline = find(contains(A,'!CAL CALIBRATION'));
    
    if ~isempty(whichline)
        CALIBRATION = split(char(A(whichline(1))));
        calibrationMethod = string(CALIBRATION{5});
        eye = string(CALIBRATION{7});
        calibrationQuality = string(CALIBRATION{8});
    
        if ~any(find(contains(A, '!CAL VALIDATION') & contains(A, 'ABORTED'))) && ~isempty(find(contains(A,'!CAL VALIDATION')))
            whichline = find(contains(A, '!CAL VALIDATION'));
            VALIDATION = split(char(A(whichline(1))));
            validationQuality = string(VALIDATION{8});
            averageError = str2double(VALIDATION{10});
            maxError = str2double(VALIDATION{12});
            offsetDeg = str2double(VALIDATION{15});
            offpix = string(VALIDATION{17});
        else
            validationQuality = "";
            averageError = NaN;
            maxError = NaN;
            offsetDeg = NaN;
            offpix = "";
        end
    else
        calibrationMethod = "";
        eye = "";
        calibrationQuality = "";
        validationQuality = "";
        averageError = NaN;
        maxError = NaN;
        offsetDeg = NaN;
        offpix = "";
    end
    
    noTrackingTrials = NaN; % fill in later
    
   % Create a struct with the data
    metaDataStruct = struct('SamplingRate', samplingRate, ...
                            'CalibrationMethod', calibrationMethod, ...
                            'Eye', eye, ...
                            'CalibrationQuality', calibrationQuality, ...
                            'ValidationQuality', validationQuality, ...
                            'AverageError', averageError, ...
                            'MaxError', maxError, ...
                            'OffsetDeg', offsetDeg, ...
                            'Offpix', offpix, ...
                            'NoTrackingTrials', noTrackingTrials, ...
                            'Filename', string(filename), ...
                            'SubjectName', string(subjectname));
    
end