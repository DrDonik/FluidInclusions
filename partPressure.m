% This routine implements Jaro's model to calculate the pressure of the
% liquid and the vapour phase. You need to run get_r first and feed
% r to this function. Feed the bubble-radius in um, I will
% convert it to SI.

function obj = partPressure(obj)

for T_ctr = find(obj.store_p_l==0)
    
    if isnan(obj.r(T_ctr))
        % It would make sense if the p_l on the iso_Th is calculated,
        % instead of giving back NaN, but water95 can't do this;
        if obj.store_T(T_ctr) <= 273.15
            obj.store_p_l(T_ctr) = NaN;
        else
            obj.store_p_l(T_ctr) = pressure(obj.store_rho_overall_at_T(T_ctr), obj.store_T(T_ctr));
            num2str(T_ctr)
        end
        obj.store_p_v(T_ctr) = NaN;
        
    else
        P0 = saturationPressure(obj.store_T(T_ctr));
        rho_liquid = inclusion.liqvap_density(obj.store_T(T_ctr));
        rho_vapour = inclusion.liqvap_density_vapour(obj.store_T(T_ctr));

        delta_P = 2*inclusion.surface_tension(obj.store_T(T_ctr))/(obj.r(T_ctr)/1e6);

        obj.store_p_l(T_ctr) = P0 - rho_liquid/(rho_liquid - rho_vapour)*delta_P;
        obj.store_p_v(T_ctr) = P0 - rho_vapour/(rho_liquid - rho_vapour)*delta_P;
    end;
end