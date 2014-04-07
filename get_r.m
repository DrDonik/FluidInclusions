% This routine tries to make use of the isochoricobjective helper function
% to find the minima of the helmholtz energy.

% dg denotes vapour-density divided by the critical density of water
% gm denotes the Bubble volume to inclusion volume ratio

% I expect T and Th_inf to be input in K, V in um (Not SI!!!)

function [r, steamDensity_corrected] = ...
    get_r(T, Th_inf, V)

% Load some data from IAPWS-95
coeffs = readIAPWS95data();

% V was entered in um^3, so change it to SI.
V = V*1e-18;

% Some constants
rhoc = 322;
Tc = 647.096;
rc = 1.1808741e-8;
b = -0.625;
mu = 1.256;

% The fit options. Usually the fit converges after less than 10 iterations,
% if it doesn't there will be no minimum. The GradObj-entry tells the fit
% routine to take the Jacobian into account.
options = optimset('TolX',1e-12,'TolFun',1e-15,'GradObj','on','Hessian','user-supplied','Algorithm','trust-region-reflective','Display','off','MaxIter',50);
rhoOverallInitial = liqvap_density(Th_inf)*1000;

% Apply the volume correction
[reftemp, alpha_V] = expansion_coeff(T);
rho_correction = ((1-(reftemp-Th_inf+273.15)*alpha_V)/(1-(reftemp-T+273.15)*alpha_V));
rho_overall_at_T = rhoOverallInitial*rho_correction;
dm = rho_overall_at_T/rhoc;

% Make an initial estimate using IAPWS-95, pretending there was no
% surface tension. These values will be larger, but close to the
% final values.
minvars_corrected(1) = (1 - rho_overall_at_T/liqvap_density(T)/1000);
minvars_corrected(2) = liqvap_density_vapor(T)/rhoc*1000;

% Calculate the surface tension
tau = Tc/T;
stprime = rc/V.^(1/3) * ((tau - 1)/tau)^mu * ((1 + b)*tau - b);

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
	r = NaN;
	if (nargout > 1)
		% The density was asked too, so save it.
		steamDensity_corrected = NaN;
	end;
	exitflag_corrected = 0;
	
end;

if ~(~exitflag_corrected || minvars_corrected(1) < 0 || isnan(minvars_corrected(1)) || ~isreal(minvars_corrected(1))  || minvars_corrected(1) <= 1e-9)

	% There is a bubble possible. Calculate it's radius and save it
	r = (3*V.*minvars_corrected(1)./ (4 * pi)).^(1/3)*1e6;
	if (nargout > 1)
		% The density was asked too, so save it.
		steamDensity_corrected = minvars_corrected(2) * rhoc;
	end;

end;
        
% Done
return;