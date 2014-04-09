% This routine implements Jaro's model to calculate the pressure of the
% liquid and the vapour phase. You need to run get_r first and feed
% r to this function. Feed the bubble-radius in um, I will
% convert it to SI.

function obj = partPressure(obj)

for T_ctr = find(obj.store_p_l==0)

    P0 = saturationPressure(obj.store_T(T_ctr));
    rho_liquid = inclusion.liqvap_density(obj.store_T(T_ctr));
    rho_vapour = inclusion.liqvap_density_vapour(obj.store_T(T_ctr));

    delta_P = 2*inclusion.surface_tension(obj.store_T(T_ctr))/(obj.r/1e6);

    obj.store_p_l = P0 - rho_liquid/(rho_liquid - rho_vapour)*delta_P;
    obj.store_p_v = P0 - rho_vapour/(rho_liquid - rho_vapour)*delta_P;

end