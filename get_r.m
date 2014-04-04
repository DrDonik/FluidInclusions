% This routine tries to make use of the isochoricobjective helper function
% to find the minima of the helmholtz energy.

% dg denotes vapour-density divided by the critical density of water
% gm denotes the Bubble volume to inclusion volume ratio

% I expect T_array and Th_inf to be input in K, V in um (Not SI!!!)

function [radius_out_corrected_helmholtz, steamDensity_corrected, ...
    radius_out_uncorrected_helmholtz, steamDensity_uncorrected] = ...
    get_r(T_array, Th_inf, V)

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

% Initialise the output arrays to zero.
radius_out_corrected_helmholtz = zeros(size(T_array));
if nargout > 1
    steamDensity_corrected = zeros(size(T_array));
    if nargout > 2
        radius_out_uncorrected_helmholtz = zeros(size(T_array));
        if nargout > 3
            steamDensity_uncorrected = zeros(size(T_array));
        end;
    end;
end;

% The fit options. Usually the fit converges after less than 10 iterations,
% if it doesn't there will be no minimum. The GradObj-entry tells the fit
% routine to take the Jacobian into account.
options = optimset('TolX',1e-12,'TolFun',1e-15,'GradObj','on','Hessian','user-supplied','Algorithm','trust-region-reflective','Display','off','MaxIter',50);
rhoOverallInitial = liqvap_density(Th_inf)*1000;

% Walk through the supplied temperatures.
for ctr = 1:length(T_array);
    T_working = T_array(ctr);

    % Apply the volume correction
    [reftemp, alpha_V] = expansion_coeff(T_working);
    rho_correction = ((1-(reftemp-Th_inf+273.15)*alpha_V)/(1-(reftemp-T_working+273.15)*alpha_V));
    rho_overall_at_T = rhoOverallInitial*rho_correction;
    dm = rho_overall_at_T/rhoc;

    % Make an initial estimate using IAPWS-95, pretending there was no
    % surface tension. These values will be larger, but close to the
    % final values.
    minvars_corrected(1) = (1 - rho_overall_at_T/liqvap_density(T_working)/1000);
    minvars_corrected(2) = liqvap_density_vapor(T_working)/rhoc*1000;
    if nargout > 2; minvars_uncorrected = minvars_corrected; end;

    % Calculate the surface tension
    tau = Tc/T_working;
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
        radius_out_corrected_helmholtz(ctr) = NaN;
        if (nargout > 1)
            % The density was asked too, so save it.
            steamDensity_corrected(ctr) = NaN;
        end;
        exitflag_corrected = 0;
        
    end;
    
    if ~(~exitflag_corrected || minvars_corrected(1) < 0 || isnan(minvars_corrected(1)) || ~isreal(minvars_corrected(1))  || minvars_corrected(1) <= 1e-9)

        % There is a bubble possible. Calculate it's radius and save it
        radius_out_corrected_helmholtz(ctr) = (3*V.*minvars_corrected(1)./ (4 * pi)).^(1/3)*1e6;
        if (nargout > 1)
            % The density was asked too, so save it.
            steamDensity_corrected(ctr) = minvars_corrected(2) * rhoc;
        end;

    end;
        
    if nargout > 2
        % The user wanted to know the bubble radius without surface tension
        stprime = 0;
        
        % This will be the funtion to minimise
        helmholtz_function = @(minvars) isochoricobjective(minvars(1), minvars(2), tau, A, stprime, coeffs, dm);
        
        % Minimise the helmholtz energy
        try
            [minvars_uncorrected, ~, exitflag] = fmincon(helmholtz_function, minvars_uncorrected, [],[],[],[],[0 0],[1 Inf],[],options);
        catch
            
            % There seems to be a problem minimising the energy. Probably,
            % there is no bubble possible
            radius_out_uncorrected_helmholtz(ctr) = NaN;
            if nargout > 3
                % The density was asked too, so save it.
                steamDensity_uncorrected(ctr) = NaN;
            end;
            exitflag = 0;

        end;

        if ~(~exitflag || minvars_uncorrected(1) < 0 || isnan(minvars_uncorrected(1)) || ~isreal(minvars_uncorrected(1))  || minvars_corrected(1) <= 1e-9)
            % There is a bubble possible. Calculate it's radius and save it
            radius_out_uncorrected_helmholtz(ctr) = (3 * V * minvars_uncorrected(1) / (4 * pi))^(1/3)*1e6;
            if nargout > 3
                % The density was asked too, so save it.
                steamDensity_uncorrected(ctr) = minvars_uncorrected(2) * rhoc;
            end;
        end;

    end;
    
end;

% For all entries that are still 0 there was no bubble possible, change
% them to NaN. (This makes plotting easier)
mask_radius_out_corrected_helmholtz = radius_out_corrected_helmholtz == 0;
radius_out_corrected_helmholtz(mask_radius_out_corrected_helmholtz) = NaN;

if nargout > 2
    mask_radius_out_uncorrected_helmholtz = radius_out_uncorrected_helmholtz == 0;
    radius_out_uncorrected_helmholtz(mask_radius_out_uncorrected_helmholtz) = NaN;
end;

% Done
return;