% This routine loads the mineral to use when doing inclusion calculations

function [mineralNumber, pressureMinimum]  = get_fi_mineral()

    global MineralNumber PressureMinimum
    
    if isempty(MineralNumber)
        
        [mineralNumber, pressureMinimum] = set_fi_mineral();
        
    else
    
        mineralNumber = MineralNumber;
        if nargout > 1; pressureMinimum = PressureMinimum; end
        
    end
    
    
return;