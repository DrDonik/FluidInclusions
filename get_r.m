% This routine tries to make use of the isochoricobjective helper function
% to find the minima of the helmholtz energy.

% dg denotes vapour-density divided by the critical density of water
% gm denotes the Bubble volume to inclusion volume ratio

% I expect T to be in K

function obj = get_r(obj)

% Load some data from IAPWS-95
coeffs = inclusion.readIAPWS95data();

% Some constants
rhoc = 322;
Tc = 647.096;
rc = 1.1808741e-8;
b = -0.625;
mu = 1.256;

% The fit options. Usually the fit converges after less than 50 iterations,
% if it doesn't there will be no minimum. The GradObj-entry tells the fit
% routine to take the Jacobian into account.
options = optimset('TolX',1e-12,'TolFun',1e-15,'GradObj','on','Hessian','user-supplied','Algorithm','trust-region-reflective','Display','off','MaxIter',50);

for T_ctr = find(obj.store_r == 0)

    % Apply the volume correction
    dm = obj.rho_overall_at_T(T_ctr)/rhoc;

    % Make an initial estimate using IAPWS-95, pretending there was no
    % surface tension. These values will be larger, but close to the
    % final values.
    minvars_corrected(1) = (1 - obj.rho_overall_at_T(T_ctr)/inclusion.liqvap_density(obj.store_T(T_ctr))/1000);
    minvars_corrected(2) = inclusion.liqvap_density_vapour(obj.store_T(T_ctr))/rhoc*1000;


    % Calculate the surface tension
    tau = Tc/obj.store_T(T_ctr);
    stprime = rc/(obj.store_V)^(1/3) * ((tau - 1)/tau)^mu * ((1 + b)*tau - b);

    % Put the salt term here: A = C*w*Mw/Ms*dm.
    % C is the dissociation number, w the weight fraction of
    % salt, Mw and Ms the molecular weight of the salt and the
    % water. This is 2*0.05*0.308*dm for 5% NaCl.
    A = 0;

    % This will be the funtion to minimise
    helmholtz_function = @(minvars) isochoricobjective(minvars(1), minvars(2), tau, A, stprime, coeffs, dm);

    % Minimise the helmholtz energy
    try
        [minvars_corrected, ~, exitflag_corrected] = fmincon(helmholtz_function, minvars_corrected, [],[],[],[],[0 0],[1 Inf],[],options);
    catch

        % There seems to be a problem minimising the energy. Probably,
        % there is no bubble possible
        obj.store_r(T_ctr) = NaN;
        return;
    end

    if ~(~exitflag_corrected || minvars_corrected(1) < 0 || isnan(minvars_corrected(1)) || ~isreal(minvars_corrected(1))  || minvars_corrected(1) <= 1e-9)

        % There is a bubble possible. Calculate it's radius and save it
        obj.store_r(T_ctr) = (3*obj.store_V.*minvars_corrected(1)./ (4 * pi)).^(1/3)*1e6;
    else
        % There seems to be a problem minimising the energy. Probably,
        % there is no bubble possible
        obj.store_r(T_ctr) = NaN;
    end
    
end

% Done
return