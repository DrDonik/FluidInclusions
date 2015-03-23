% This routine implements Jaro's model to calculate the pressure of the
% liquid and the vapour phase. If there is no bubble possible, the function
% returns the pressure on the iso-Th-curve

function partPressure(obj)

    for T_ctr = find(obj.store_p_l==0)

        if isnan(obj.r(T_ctr))
            obj.store_p_l(T_ctr) = directPressureRaw(obj.store_rho_overall_at_T(T_ctr), obj.store_T(T_ctr));
            obj.store_p_isoTh(T_ctr) = obj.store_p_l(T_ctr);
        else
            P0 = saturationPressure(obj.store_T(T_ctr));
            rho_liquid = inclusion.liqvap_density(obj.store_T(T_ctr));
            rho_vapour = inclusion.liqvap_density_vapour(obj.store_T(T_ctr));

            delta_P = 2*inclusion.surface_tension(obj.store_T(T_ctr))/(obj.r(T_ctr)/1e6);

            obj.store_p_l(T_ctr) = P0 - rho_liquid/(rho_liquid - rho_vapour)*delta_P;
            obj.store_p_v(T_ctr) = P0 - rho_vapour/(rho_liquid - rho_vapour)*delta_P;
            obj.store_p_isoTh(T_ctr) = directPressureRaw(obj.store_rho_overall_at_T(T_ctr), obj.store_T(T_ctr));
        end
    end

    return
end
