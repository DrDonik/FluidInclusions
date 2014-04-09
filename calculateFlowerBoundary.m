% Th_obs_is_T_bin has to be set either to true (or 1, the default) for bin or to false
% (or 0) for sp.

function [Th_inf_working, r] = calculateFlowerBoundary(obj, Th_obs_is_T_bin)

% Load some data from IAPWS-95
coeffs = inclusion.readIAPWS95data();

if nargin == 1; Th_obs_is_T_bin = 1; end;

% The fit options. Usually the fit converges after less than 10 iterations,
% if it doesn't there will be no minimum. The GradObj-entry tells the fit
% routine to take the Jacobian into account.
TolX = 1e-13;
TolFun = 1e-15;
options = optimset('TolX',TolX,'TolFun',TolFun,'GradObj','on','Hessian','user-supplied','Algorithm','trust-region-reflective','Display','off','MaxIter',10);

% Some constants
rhoc = 322;
Tc = 647.096;
rc = 1.1808741e-8;
b = -0.625;
mu = 1.256;

step = 1.25;
Th_inf_working = 20 + 273.15;

while step >= 0.002/5^5

    gm_out_corrected = 1;

    while gm_out_corrected > 1e-9

        Th_inf_working = Th_inf_working - step;
        while Th_inf_working <= obj.store_T_pressureMinimum
            step = step/5;
            Th_inf_working = obj.store_T_pressureMinimum + step;
        end;

        rhoOverallInitial = inclusion.liqvap_density(Th_inf_working)*1000;

        % Apply the volume correction
        [reftemp, alpha_V] = expansion_coeff(obj, Th_inf_working);
        rho_overall_at_T = rhoOverallInitial*((1-(reftemp-Th_inf_working+273.15)*alpha_V)/(1-(reftemp-obj.store_T_pressureMinimum+273.15)*alpha_V));
        dm = rho_overall_at_T/rhoc;

        % Calculate the surface tension
        tau = Tc/obj.store_T_pressureMinimum;
        stprime = rc/(obj.store_V)^(1/3) * ((tau - 1)/tau)^mu * ((1 + b)*tau - b);

        % Put the salt term here: A = C*w*Mw/Ms*dm.
        % C is the dissociation number, w the weight fraction of
        % salt, Mw and Ms the molecular weight of the salt and the
        % water. This is 2*0.05*0.308*dm for 5% NaCl.
        A = 0;

        % This will be the funtion to minimise
        helmholtz_function = @(minvars) isochoricobjective(minvars(1), minvars(2), tau, A, stprime, coeffs, dm);

        % Make an initial estimate using IAPWS-95, pretending there was no
        % surface tension. These values will be larger, but close to the
        % final values.
        minvars_corrected(1) = (1 - rho_overall_at_T/1000/inclusion.liqvap_density(obj.store_T_pressureMinimum));
        minvars_corrected(2) = inclusion.liqvap_density_vapour(obj.store_T_pressureMinimum)/rhoc*1000;

        if minvars_corrected(1) > 0
            % Minimise the helmholtz energy
            minvars_corrected = fmincon(helmholtz_function, minvars_corrected, [],[],[],[],[0 0],[1 Inf],[],options);
        else
            continue;
        end;

        if Th_obs_is_T_bin && helmholtz_function(minvars_corrected) > helmholtz_function([0, minvars_corrected(2)])
            gm_out_corrected = 0;
        else
            gm_out_corrected = minvars_corrected(1);
        end;

    end;

    Th_inf_working = Th_inf_working + step;
    step = step/5;

end

r = (3 * obj.store_V * gm_out_corrected / (4 * pi))^(1/3)*1e6;

return
