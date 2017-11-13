% This routine implements Jaro's model to calculate the pressure of the
% liquid and the vapour phase. If there is no bubble possible, the function
% returns the pressure on the iso-Th-curve

function partPressure(obj)

    for T_ctr = find(obj.store_p_l==0)

        if isnan(obj.r(T_ctr))
            obj.store_p_l(T_ctr) = directPressureRaw(obj.store_rho_overall_at_T(T_ctr), obj.store_T(T_ctr));
            obj.store_p_v(T_ctr) = NaN;
            obj.store_p_isoTh(T_ctr) = obj.store_p_l(T_ctr);
            obj.store_p_s(T_ctr) = saturationPressure(obj.store_T(T_ctr));
        else
            [P0, liqvap_liqrho, liqvap_vaprho] = saturationPressure(obj.store_T(T_ctr));
            delta_P = 2*inclusion.surface_tension(obj.store_T(T_ctr))/(obj.r(T_ctr)/1e6);

            obj.store_p_l(T_ctr) = P0 - liqvap_liqrho/(liqvap_liqrho - liqvap_vaprho)*delta_P;
            obj.store_p_v(T_ctr) = P0 - liqvap_vaprho/(liqvap_liqrho - liqvap_vaprho)*delta_P;
            obj.store_p_isoTh(T_ctr) = directPressureRaw(obj.store_rho_overall_at_T(T_ctr), obj.store_T(T_ctr));
            obj.store_p_s(T_ctr) = P0;
        end
    end

    return
end
