% This routine sets the mineral to use when doing inclusion calculations

function [mineralNumber, T_pressureMinimum, mineral] = set_fi_mineral(mineralNumber)
    
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
            if isempty(mineralNumber); mineralNumber = 7; end
        else
            mineral_list = {'Calcite, reference T = 30°, Rao', ...
                'Calcite, reference T = 20° Rao, Maximum Density at 5.15 °C', ...
                'Calcite, reference T = 20°, TPPM, Maximum Density at 4.65 °C', ...
                'Quartz, reference T = 20° TPPM, Maximum Density at 6.05 °C', ...
                'Quartz, reference T = 20° HTMIAC, Maximum Density at 6.15 °C', ...
                'Gypsum, reference T = 19.4° Schofield, Maximum Density at 8.35 °C', ...
                'No correction'};
            [mineralNumber,isok] = listdlg('PromptString','Select a Mineral (invoke set_fi_mineral to change it later):','SelectionMode','single','ListSize',[500,200],'ListString',mineral_list);
            if ~isok; mineralNumber = 7; end
        end
    end
    
    switch mineralNumber
        case 1; T_pressureMinimum = 5.15 + 273.15; mineral = 'Calcite';
        case 2; T_pressureMinimum = 5.15 + 273.15; mineral = 'Calcite (Rao)';
        case 3; T_pressureMinimum = 4.65 + 273.15; mineral = 'Calcite (TPPM)';
        case 4; T_pressureMinimum = 6.05 + 273.15; mineral = 'Quartz (TPPM)';
        case 5; T_pressureMinimum = 6.15 + 273.15; mineral = 'Quartz (HTMIAC)';
        case 6; T_pressureMinimum = 8.35 + 273.15; mineral = 'Gypsum (Schofield)';
        case 7; T_pressureMinimum = 4.00 + 273.15; mineral = 'No correction';
    end

return