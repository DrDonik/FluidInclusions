% Th_obs_is_T_bin has to be set either to true (or 1, the default) for bin or to false
% (or 0) for sp.

function [Th_inf, bubble_radius_out] = flower_boundary(V, Th_obs_is_T_bin)

% Check for the host mineral. If you want to change it, run set_fi_mineral
[~, pressureMinimum] = get_fi_mineral();

% Load some data from IAPWS-95
coeffs = readIAPWS95data();

if nargin == 1; Th_obs_is_T_bin = 1; end;

% The fit options. Usually the fit converges after less than 10 iterations,
% if it doesn't there will be no minimum. The GradObj-entry tells the fit
% routine to take the Jacobian into account.
TolX = 1e-13;
TolFun = 1e-15;
options = optimset('TolX',TolX,'TolFun',TolFun,'GradObj','on','Hessian','user-supplied','Algorithm','trust-region-reflective','Display','off','MaxIter',50);

% Some constants
rhoc = 322;
Tc = 647.096;
rc = 1.1808741e-8;
b = -0.625;
mu = 1.256;

% V was entered in um³, so change it to SI.
V = V*1e-18;
Th_inf = (70+273.15)*ones(size(V));
gm_out_corrected = zeros(size(V));

for ctr = 1:length(V)

    step = 1.25;

    while step >= 0.002/5^5

        gm_out_corrected(ctr) = 1;

        while gm_out_corrected(ctr) > 1e-9

            Th_inf(ctr) = Th_inf(ctr) - step;
            while Th_inf(ctr) <= pressureMinimum
                step = step/5;
                Th_inf(ctr) = pressureMinimum+step;
            end;

            rhoOverallInitial = liqvap_density(Th_inf(ctr))*1000;

            % Apply the volume correction
            [reftemp, alpha_V] = expansion_coeff(Th_inf(ctr));
            rho_overall_at_T = rhoOverallInitial*((1-(reftemp-Th_inf(ctr)+273.15)*alpha_V)/(1-(reftemp-pressureMinimum+273.15)*alpha_V));
            dm = rho_overall_at_T/rhoc;

            % Calculate the surface tension
            tau = Tc/pressureMinimum;
            stprime = rc/V^(1/3) * ((tau - 1)/tau)^mu * ((1 + b)*tau - b);

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
            minvars_corrected(1) = (1 - rho_overall_at_T/1000/liqvap_density(pressureMinimum));
            minvars_corrected(2) = liqvap_density_vapor(pressureMinimum)/rhoc*1000;

            if minvars_corrected(1) > 0
                % Minimise the helmholtz energy
                minvars_corrected = fmincon(helmholtz_function, minvars_corrected, [],[],[],[],[0 0],[1 Inf],[],options);
            else
                continue;
            end;

            if Th_obs_is_T_bin && helmholtz_function(minvars_corrected) > helmholtz_function([0, minvars_corrected(2)])
                gm_out_corrected(ctr) = 0;
            else
                gm_out_corrected(ctr) = minvars_corrected(1);
            end;

        end;

        Th_inf(ctr) = Th_inf(ctr) + step;
        step = step/5;

    end
    
end

bubble_radius_out = (3 * V * gm_out_corrected / (4 * pi))^(1/3)*1e6;

return
