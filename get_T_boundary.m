%get_T_boundary  Calculate the boundary temperature for inclusions with
%specific volumes and densities
%
% [T_boundary, r_boundary] = 
%     get_T_boundary(Th_inf, V, calc_sp_boundary, calc_prograde_boundary, pressureMinimum)
%
% calculates the boundary temperature of an inclusion of a given size and a
% given nominal homogenisation temperature.
%   Th_inf is the nominal homogenization temperature (in Kelvin)
%   V is the volume (in micrometers^3) for which the boundary temperatures 
%     should be calculated
%   calc_sp_boundary is whether you want the spinodal (1) or the binodal (0) 
%     temperature
%   calc_prograde_boundary is whether you want the prograde (1) or retrograde (0) 
%     temperature
%   pressureMinimum is the temperature of minimum pressure in the inclusion
%   (and is used as is the initial guess for T_boundary)
% All arguments are mandatory.
% 
% The output contains the calculated boundary temperature, T_boundary (in K) 
% and the calculated bubble radius at this temperatures, r_boundary (in um).
%

function [T_boundary, r_boundary] = ...
    get_T_boundary(Th_inf, V, calc_sp_boundary, calc_prograde_boundary, pressureMinimum)

% Load some data from IAPWS-95
coeffs = inclusion.readIAPWS95data();

tolerance = inclusion.get_tolerance;

% The fit options. Usually the fit converges after less than 10 iterations,
% if it doesn't there will be no minimum. The GradObj-entry tells the fit
% routine to take the Jacobian into account.
TolX = 1e-13;
TolFun = 1e-15;
options = optimset('TolX',TolX,'TolFun',TolFun,'GradObj','on',...
    'Hessian','user-supplied','Algorithm','trust-region-reflective',...
    'Display','off','MaxIter',10);

% Some constants
rhoc = 322;
Tc = 647.0960;
rc = 1.1808741e-8;
b = -0.625;
mu = 1.256;

% V was entered in um^3, so change it to SI.
V = V*1e-18;

rhoOverallInitial = inclusion.liqvap_density(Th_inf)*1000;

if calc_prograde_boundary; dir = 1; else dir = -1; end;
T_boundary_working = pressureMinimum;

step = 1.25;

iterationCounter = 0;

while step >= tolerance
    % We stop as soon as we reach the tolerance in temperature.
    
    iterationCounter = iterationCounter + 1;
    
    if iterationCounter > 1
        % We crossed T_boundary; make the step smaller and change direction.
        step = step/5;
        dir = -dir;
    end;
    
    while iterationCounter == 1 || isnan(gm_out_corrected) == isnan(previous_gm)
        % As long as we don't cross T_boundary, the above condition
        % will hold. We keep on heating (or cooling for retrograde)
        % the inclusion and check if there is still a bubble
        
        if (sign(T_boundary_working + dir*step - pressureMinimum) ~= sign(T_boundary_working - pressureMinimum)) ...
                && ((T_boundary_working - pressureMinimum) ~= 0)
            % We will cross the pressure minimum. This is bad. So
            % change to the pressure minimum.
            if isnan(gm_out_corrected)
                T_boundary_working = pressureMinimum;
                disp('Returned to the pressure minimum')
                gm_out_corrected = 1;
                break;
            else
                disp('Sending Oops! I wanted to go to the pressure minimum from a non-NaN-radius')
                keyboard
            end
        else
            while sign(T_boundary_working + dir*step - Th_inf) == sign(calc_prograde_boundary);
                % Make sure we don't cross Th_inf, since we know we
                % will not find a bubble beyond that temperature.
                step = step/5;
            end;
            T_boundary_working = T_boundary_working + dir*step;
        end;
        
        % Apply the volume correction
        [reftemp, alpha_V] = inclusion.expansion_coeff(T_boundary_working);
        rho_overall_at_T = rhoOverallInitial*((1-(reftemp-Th_inf+273.15)*alpha_V)/(1-(reftemp-T_boundary_working+273.15)*alpha_V));
        dm = rho_overall_at_T/rhoc;
        
        % Calculate the surface tension
        tau = Tc/T_boundary_working;
        stprime = rc/V^(1/3) * ((tau - 1)/tau)^mu * ((1 + b)*tau - b);
        
        % Put the salt term here: A = C*w*Mw/Ms*dm.
        % C is the dissociation number, w the weight fraction of
        % salt, Mw and Ms the molecular weight of the salt and the
        % water. This is 2*0.05*0.308*dm for 5% NaCl.
        A = 0;
        
        % This will be the funtion to minimise
        helmholtz_function = @(minvars) inclusion.isochoricobjective(minvars(1), minvars(2), tau, A, stprime, coeffs, dm);
        
        % Make an initial estimate using IAPWS-95, pretending there was no
        % surface tension. These values will be larger, but close to the
        % final values.
        minvars_corrected(1) = (1 - rho_overall_at_T/1000/inclusion.liqvap_density(T_boundary_working));
        minvars_corrected(2) = inclusion.liqvap_density_vapour(T_boundary_working)/rhoc*1000;
        
        if minvars_corrected(1) > 0
            % Minimise the helmholtz energy
            minvars_corrected = fmincon(helmholtz_function, minvars_corrected, [],[],[],[],[0 0],[1 Inf],[],options);
            gm_out_corrected = minvars_corrected(1);
            
            if ~isreal(gm_out_corrected) || gm_out_corrected <= 1e-9
                % There was no bubble possible at this temperature and this volume.
                gm_out_corrected = NaN;
            elseif ~calc_sp_boundary && helmholtz_function(minvars_corrected) > helmholtz_function([0, minvars_corrected(2)])
                % The minimum was found, but it is already in the metastable region
                gm_out_corrected = NaN;
            end;
        else
            gm_out_corrected = NaN;
        end;
        
        if iterationCounter == 1; break; end;
        
    end;
    % When we reach here, we crossed T_boundary.
    
    if iterationCounter == 1;
        if isnan(gm_out_corrected)
            % No minimum found; most probably there is no bubble possible
            % in this volume at this temperature;
            T_boundary = NaN;
            r_boundary = NaN;
            
            return;
        end;
    else
        temp_save_gm = previous_gm;
    end;
    
    % Save the radius for the next iteration.
    previous_gm = gm_out_corrected;
    
end;

% We found a T_boundary (or gave up). Right now we have to find out
% whether we're on the NaN-side of T_boundary or on the side with a
% bubble, and then we'll save the value for T_boundary accordingly.
if isnan(gm_out_corrected)
    T_boundary = T_boundary_working-dir*step;
    r_boundary = (3 * V * temp_save_gm / (4 * pi))^(1/3)*1e6;
else
    T_boundary = T_boundary_working;
    r_boundary = (3 * V * gm_out_corrected / (4 * pi))^(1/3)*1e6;
end;

% Done.
return;
