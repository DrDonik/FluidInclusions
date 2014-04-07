%get_T_bin  Calculate the binodal temperature for inclusions with
%specific volumes and densities
%
% [T_bin, r_bin] = 
%                   get_T_bin(Th_inf, V, T_bin_to_calc, pressureMinimum)
%
% calculates the binodal temperature of an inclusion of a given size and a
% given nominal homogenisation temperature.
% 	Th_inf is the nominal homogenization temperature (in Kelvin)
%   V is the volume (in micrometers^3) for which the binodal temperatures 
%     should be calculated
%   T_sp_to_calc is whether the user wants prograde (1) or retrograde (-1) 
%     temperatures
%   pressureMinimum is the temperature of minimum pressure in the inclusion
%   (and is used as is the initial guess for T_sp)
% All arguments are mandatory.
% 
% The output contains the calculated binodal temperature, T_bin (in K) 
% and the calculated bubble radius at this temperatures, r_bin (in um).
%

function [T_bin, r_bin] = ...
    get_T_bin(Th_inf, V, T_bin_to_calc, pressureMinimum)

% Load some data from IAPWS-95
coeffs = readIAPWS95data();

tolerance = get_tolerance;

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

% V was entered in um^3, so change it to SI.
V = V*1e-18;

rhoOverallInitial = liqvap_density(Th_inf)*1000;

dir = T_bin_to_calc;
T_bin_working = pressureMinimum;

step = 1.25;

iterationCounter = 0;

while step/5 >= tolerance
	% We stop as soon as we reach the tolerance in temperature.
	
	iterationCounter = iterationCounter + 1;
	
    if iterationCounter > 1
		% We crossed T_bin; make the step smaller and change direction.
		step = step/5;
		dir = -dir;
	end;
	
	while iterationCounter == 1 || isnan(gm_out_corrected) == isnan(previous_gm)
		% As long as we don't cross T_bin, the above condition
		% will hold. We keep on heating (or cooling for retrograde)
		% the inclusion and check if there is still a stable bubble

		if (sign(T_bin_working + dir*step - pressureMinimum) ~= sign(T_bin_working - pressureMinimum)) ...
				&& ((T_bin_working - pressureMinimum) ~= 0)
            % We will cross the pressure minimum. This is bad. So
            % change to the pressure minimum.
            if isnan(gm_out_corrected)
                T_bin_working = pressureMinimum;
                disp('Returned to the pressure minimum')
                gm_out_corrected = 1;
                break;
            else
                disp('Sending Oops! I wanted to go to the pressure minimum from a non-NaN-radius')
                keyboard
            end
		else
            while sign(T_bin_working + dir*step - Th_inf) == sign(T_bin_to_calc);
                % Make sure we don't cross Th_inf, since we know we
                % will not find a bubble beyond that temperature.
                step = step/5;
            end;
			T_bin_working = T_bin_working + dir*step;
		end;
		
		% Apply the volume correction
		[reftemp, alpha_V] = expansion_coeff(T_bin_working);
		rho_overall_at_T = rhoOverallInitial*((1-(reftemp-Th_inf+273.15)*alpha_V)/(1-(reftemp-T_bin_working+273.15)*alpha_V));
		dm = rho_overall_at_T/rhoc;
		
		% Calculate the surface tension
		tau = Tc/T_bin_working;
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
		minvars_corrected(1) = (1 - rho_overall_at_T/1000/liqvap_density(T_bin_working));
		minvars_corrected(2) = liqvap_density_vapor(T_bin_working)/rhoc*1000;
		
		if minvars_corrected(1) > 0
			% Minimise the helmholtz energy
			minvars_corrected = fmincon(helmholtz_function, minvars_corrected, [],[],[],[],[0 0],[1 Inf],[],options);
			gm_out_corrected = minvars_corrected(1);
		
			if ~isreal(gm_out_corrected) || gm_out_corrected <= 1e-9
				% There was no bubble possible at this temperature int
				% his volume.
				gm_out_corrected = NaN;
			elseif helmholtz_function(minvars_corrected) > helmholtz_function([0, minvars_corrected(2)])
				% The minimum was found, but it is already in the metastable
				% region
				gm_out_corrected = NaN;
			end;
		else
			gm_out_corrected = NaN;
		end;
		
        if iterationCounter == 1; break; end;
        
	end;
	% When we reach here, we crossed T_bin.
	
    if iterationCounter == 1;
        if isnan(gm_out_corrected)
            % No minimum found; most probably there is no bubble possible
            % in this volume at this temperature;
            T_bin = NaN;
            r_bin = NaN;

            return;
        end;
    else
        temp_save_gm = previous_gm;
    end;

    % Save the radius for the next iteration.
    previous_gm = gm_out_corrected;

end;

% We found a T_bin (or gave up). Right now we have to find out
% whether we're on the NaN-side of T_bin or on the side with a
% stable bubble, and then we'll save the value for T_bin accordingly.
if isnan(gm_out_corrected)
	T_bin = T_bin_working-dir*step;
	r_bin = (3 * V * temp_save_gm / (4 * pi))^(1/3)*1e6;
else
	T_bin = T_bin_working;
	r_bin = (3 * V * gm_out_corrected / (4 * pi))^(1/3)*1e6;
end;

% Done.
return;
