%get_T_sp  Calculate the spinodal temperature for inclusions with
%specific volumes and densities (multicore version)
%
% [T_sp, r_sp] = 
%              get_T_sp_multicore(Th_inf, V, T_sp_to_calc, pressureMinimum)
%
% calculates the spinodal temperature of an inclusion of a given size and
% nominal homogenisation temperature.
%   Th_inf is the nominal homogenization temperature (in Kelvin)
%   V is the volume (in micrometers^3) for which the spinodal temperatures 
%     should be calculated
%   T_sp_to_calc is whether the user wants prograde (1) or retrograde (-1) 
%     temperatures
%   pressureMinimum is the temperature of minimum pressure in the inclusion
%   (and is used as is the initial guess for T_sp)
% All arguments are mandatory.
% 
% The output contains the calculated spinodal temperature, T_sp (in K) 
% and the calculated bubble radius at this temperatures, r_sp (in um).
%

function [T_sp, r_sp] = ...
    get_T_sp(Th_inf, V, T_sp_to_calc, pressureMinimum)

% Load some data from IAPWS-95
coeffs = readIAPWS95data();

% The fit options. Usually the fit converges after less than 10 iterations,
% if it doesn't there will be no minimum. The GradObj-entry tells the fit
% routine to take the Jacobian into account.
TolX = 1e-13;
TolFun = 1e-15;
options = optimset('TolX',TolX,'TolFun',TolFun,'GradObj','on',...
    'Hessian','user-supplied','Algorithm','trust-region-reflective',...
    'Display','off','MaxIter',10);

% Some constants
rhoc = coeffs.rhoc;
Tc = coeffs.Tc;
rc = 1.1808741e-8;
b = -0.625;
mu = 1.256;

rhoOverallInitial = liqvap_density(Th_inf)*1000;

dir = T_sp_to_calc;
T_sp_working = pressureMinimum;

step = 1.25;

iterationCounter = 0;

while step > tolerance
    % We stop as soon as we reached the tolerance in temperature.
    
    iterationCounter = iterationCounter + 1;

    if iterationCounter > 1
        % We crossed T_sp; make the step smaller and change direction.
        step = step/5;
        dir = -dir;
    end;
    
    while iterationCounter == 1 || isnan(radius_out_corrected) == isnan(previous_radius)
        % As long as we don't cross T_sp, the above condition
        % will hold. We keep on heating (or cooling for retrograde)
        % the inclusion and check if there is still a bubble
        
        if (sign(T_sp_working + dir*step - pressureMinimum) ~= sign(T_sp_working - pressureMinimum)) ...
                && ((T_sp_working - pressureMinimum) ~= 0)
            % We will cross the pressure minimum. This is bad. So
            % change to the pressure minimum.
            if isnan(radius_out_corrected)
                T_sp_working = pressureMinimum;
                disp('Returned to the pressure minimum')
                radius_out_corrected = 1;
                break;
            else
                disp('Sending Oops! I wanted to go to the pressure minimum from a non-NaN-radius')
                keyboard
            end
        else
            while sign(T_sp_working + dir*step - Th_inf) == sign(T_sp_to_calc);
                % Make sure we don't cross Th_inf, since we know we
                % will not find a bubble beyond that temperature.
                step = step/5;
            end;
            T_sp_working = T_sp_working + dir*step;
        end;

        % Apply the volume correction
        [reftemp, alpha_V] = expansion_coeff(T_sp_working);
        rho_overall_at_T = rhoOverallInitial*((1-(reftemp-Th_inf+273.15)*alpha_V)/(1-(reftemp-T_sp_working+273.15)*alpha_V));
        dm = rho_overall_at_T/rhoc;

        % Calculate the surface tension
        tau = Tc/T_sp_working;
        stprime = rc/V^(1/3) * (((tau - 1)/tau)^mu * ((1 + b)*tau - b));

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
        minvars_corrected(1) = (1 - rho_overall_at_T/1000/liqvap_density(T_sp_working));
        minvars_corrected(2) = liqvap_density_vapor(T_sp_working)/rhoc*1000;

        if minvars_corrected(1) > 0
            % Minimise the helmholtz energy
            minvars_corrected = fmincon(helmholtz_function, minvars_corrected, [],[],[],[],[0 0],[1 Inf],[],options);
            radius_out_corrected = minvars_corrected(1);
            if ~isreal(radius_out_corrected) || radius_out_corrected <= 1e-9
                % There was no bubble possible at this temperature int
                % his volume.
                radius_out_corrected = NaN;
            end;
        else
            radius_out_corrected = NaN;
        end
        
        if iterationCounter == 1; break; end;
        
    end;
    % When we reach here, we crossed T_sp.

    if iterationCounter == 1;
        if isnan(radius_out_corrected)
            % No minimum found; most probably there is no bubble possible
            % in this volume at this temperature;
            T_sp = NaN;
            r_sp = NaN;

            return;
        end;
    else
        temp_save_radius = previous_radius;
    end;

    % Save the radius for the next iteration.
    previous_radius = radius_out_corrected;

end;

% We found a T_sp (or gave up). Right now we have to find out
% whether we're on the NaN-side of T_sp or on the side with a
% bubble, and then we'll save the value for T_sp accordingly.
if isnan(radius_out_corrected)
    T_sp = T_sp_working-dir*step;
    r_sp = (3 * V * temp_save_radius / (4 * pi))^(1/3)*1e6;
else
    T_sp = T_sp_working;
    r_sp = (3 * V * radius_out_corrected / (4 * pi))^(1/3)*1e6;
end;

% Done.
return;
