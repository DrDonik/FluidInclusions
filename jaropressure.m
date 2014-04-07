% This routine implements Jaro's model to calculate the pressure of the
% liquid and the vapour phase. You need to run get_r first and feed
% r to this function. Feed the bubble-radius in um, I will
% convert it to SI.

function [P_vapour, P_liquid] = jaropressure(T, r)

P0 = liqvap(T);
rho_liquid = liqvap_density(T);
rho_vapor = liqvap_density_vapor(T);

delta_P = 2*surface_tension(T)/(r/1e6);

P_liquid = P0-rho_liquid/(rho_liquid-rho_vapor)*delta_P;
P_vapour = P0-rho_vapor/(rho_liquid-rho_vapor)*delta_P;

