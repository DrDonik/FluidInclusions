function tolerance = get_tolerance()
    
    global ToleranceInStore
    
    if ~isempty(ToleranceInStore)
        
        tolerance = ToleranceInStore;

    else

        if ~isempty(strfind(computer,'GLNX')) || ~isempty(strfind(computer,'linux')) || ~isempty(strfind(computer,'MAC'))
            % Linux OS or Mac OS X seems to be our environment
            PATH = '/tmp';
            USER = getenv('USER');
        elseif ~isempty(strfind(computer,'PCWIN'))
            % Windows seems to be our environment
            PATH = getenv('TEMP');
            USER = getenv('USERNAME');
        end;

        if ~exist([PATH, '/', USER, '_tolerance.mat'],'file')
            tolerance = set_tolerance();
        else
            load([PATH, '/', USER, '_tolerance.mat']);
        end;
        
    end;

return;