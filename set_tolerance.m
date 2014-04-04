function tolerance = set_tolerance(tolerance)
    
    global ToleranceInStore
    
    if ~isempty(strfind(computer,'GLNX')) || ~isempty(strfind(computer,'linux')) || ~isempty(strfind(computer,'MAC'))
        % Linux OS or Mac OS X seems to be our environment
        PATH = '/tmp';
        USER = getenv('USER');
    elseif ~isempty(strfind(computer,'PCWIN'))
        % Windows seems to be our environment
        PATH = getenv('TEMP');
        USER = getenv('USERNAME');
    end;

    if exist('tolerance','var');
        
        save([PATH, '/', USER, '_tolerance.mat'],'tolerance')
        
    else

        while ~exist([PATH, '/', USER, '_tolerance.mat'],'file')
            tolerance = input('Please indicate the tolerance to be used [1e-5]: ');
            if isempty(tolerance); tolerance = 1e-5; end;
            save([PATH, '/', USER, '_tolerance.mat'],'tolerance')
        end;
    
    end;
    
    ToleranceInStore = tolerance;

return;