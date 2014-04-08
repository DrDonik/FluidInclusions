% This routine implements Jaro's model to calculate the pressure of the
% liquid and the vapour phase. You need to run get_r first and feed
% r to this function. Feed the bubble-radius in um, I will
% convert it to SI.

function [P_vapour, P_liquid] = partPressure(obj)

P0 = saturationPressure(obj.T);
rho_liquid = inclusion.liqvap_density(obj.T);
rho_vapour = inclusion.liqvap_density_vapour(obj.T);

delta_P = 2*inclusion.surface_tension(obj.T)/(obj.r/1e6);

P_liquid = P0 - rho_liquid/(rho_liquid - rho_vapour)*delta_P;
P_vapour = P0 - rho_vapour/(rho_liquid - rho_vapour)*delta_P;

