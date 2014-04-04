% This routine implements Jaro's model to calculate the pressure of the
% liquid and the vapour phase. You need to run get_r first and feed
% bubbleRadius to this function. Feed the bubble-radius in um, I will
% convert it to SI.

function [P_vapor,P_liquid] = jaropressure(T_array, bubbleRadius)

P_liquid = zeros(size(T_array));
P_vapor = zeros(size(T_array));

ctr = 0;
for T_working = T_array
    ctr = ctr + 1;
    
    P0 = liqvap(T_working);
    rho_liquid = liqvap_density(T_working);
    rho_vapor = liqvap_density_vapor(T_working);
    
    delta_P = 2*surface_tension(T_working)/(bubbleRadius(ctr)/1e6);
    
    P_liquid(ctr) = P0-rho_liquid/(rho_liquid-rho_vapor)*delta_P;
    P_vapor(ctr) = P0-rho_vapor/(rho_liquid-rho_vapor)*delta_P;
        
end;