% This routine sets the mineral to use when doing inclusion calculations

function [mineralNumber, pressureMinimum] = set_fi_mineral(mineralNumber)

    global MineralNumber PressureMinimum;
    
    if nargin == 0 || isempty(mineralNumber)
        if ~usejava('desktop');
            mineral_list = ['1: Calcite, reference T = 30°, Rao\n', ...
                '2: Calcite, reference T = 20° Rao, Maximum Density at 5.15 °C\n', ...
                '3: Calcite, reference T = 20°, TPPM, Maximum Density at 4.65 °C\n', ...
                '4: Quartz, reference T = 20° TPPM, Maximum Density at 6.05 °C\n', ...
                '5: Quartz, reference T = 20° HTMIAC, Maximum Density at 6.15 °C\n', ...
                '6: Gypsum, reference T = 19.4° Schofield, Maximum Density at 8.35 °C\n', ...
                '7: No correction'];
            mineralNumber = input(['Select a Mineral (invoke set_fi_mineral to change it later):\n', mineral_list, ...
                '\nEnter the number corresponding to the mineral you want [7]: ']);
            if isempty(mineralNumber); mineralNumber=7; end
        else
            mineral_list = {'Calcite, reference T = 30°, Rao', ...
                'Calcite, reference T = 20° Rao, Maximum Density at 5.15 °C', ...
                'Calcite, reference T = 20°, TPPM, Maximum Density at 4.65 °C', ...
                'Quartz, reference T = 20° TPPM, Maximum Density at 6.05 °C', ...
                'Quartz, reference T = 20° HTMIAC, Maximum Density at 6.15 °C', ...
                'Gypsum, reference T = 19.4° Schofield, Maximum Density at 8.35 °C', ...
                'No correction'};
            [mineralNumber,isok] = listdlg('PromptString','Select a Mineral (invoke set_fi_mineral to change it later):','SelectionMode','single','ListSize',[500,200],'ListString',mineral_list);
            if ~isok; mineralNumber=7; end
        end
    end
    
    % Checking what OS we're running on ...
    if ~isempty(strfind(computer,'GLNX'))  || ~isempty(strfind(computer,'linux')) || ~isempty(strfind(computer,'MAC'))
        % Linux OS seems to be our environment
        PATH = '/tmp';
        USER = getenv('USER');
    elseif strfind(computer,'PCWIN')
        % Windows seems to be our environment
        PATH = getenv('TEMP');
        USER = getenv('USERNAME');
    end
    
    save([PATH, '/', USER, '_inclusion_mineral.mat'],'mineralNumber')
        
    MineralNumber = mineralNumber;
    
    switch mineralNumber
        case 1; pressureMinimum = 5.15 + 273.15; % Calcite
        case 2; pressureMinimum = 5.15 + 273.15; % Calcite Rao
        case 3; pressureMinimum = 4.65 + 273.15; % Calcite TPPM
        case 4; pressureMinimum = 6.05 + 273.15; % Quartz TPPM
        case 5; pressureMinimum = 6.15 + 273.15; % Quartz HTMIAC
        case 6; pressureMinimum = 8.35 + 273.15; % Gypsum Schofield
        case 7; pressureMinimum = 4.00 + 273.15; % No correction
    end
    
    PressureMinimum = pressureMinimum;

return