%get_T_sp_multicore  Calculate the spinodal temperature for inclusions with
%specific volumes and densities (multicore version)
%
% OutputCell = 
%              get_T_sp_multicore(InputCell)
%
% calculates the spinodal temperature of an inclusion of a given size and
% nominal homogenisation temperature. InputCell consists of 4 entries:
%   InputCell(1): Th_inf is an array of nominal homogenization temperatures
%       (in Kelvin)
%   InputCell(2): V is an array of volumes (in micrometers³) for which the
%       spinodal temperatures should be calculated
%   InputCell(3): T_sp_to_calc is whether the user wants prograde (1) or 
%       retrograde (-1) temperatures.
%   InputCell(4): The initial guess for T_sp
% All arguments are mandatory.
% 
% You will get one tupel [Th_sp(i), r_sp(i)] per cell entry (OutputCell{i}).
% The tupels contain the calculated spinodal temperatures in Th_sp (in K) 
% and the calculated bubble radii at this temperatures in r_sp (in um).
%
% This is the multicore version of <a href="matlab: help
% get_T_sp">get_T_sp</a>. To use it, prepare your input cell to contain
% one [Th_inf, V, tolerance, T_sp_to_calc] tupel per cell entry:
%     InputCell{i} = [Th_inf(i), V(i), T_sp_to_calc(i), initial_guess_T_sp(i)];
%
% Then start as many Matlab Sessions as you have processor cores and run in
% all but one of them
%     startmulticoreslave
% Now start in the last Matlab session
%     OutputCell = startmulticoremaster(@get_T_sp_multicore, InputCell);
%
% If you want to run the program on multiple machines, run in all but one
% of the Matlab sessions
%     startmulticoreslave('/path/to/a/commonly/accessible/directory/')
% Then enter in the last Matlab session
%     settings.multicoreDir = '/path/to/a/commonly/accessible/directory/';
% and start the calculation with
%     OutputCell = startmulticoremaster(@get_T_sp_multicore, InputCell, settings);
%
% See also <a href="matlab: help get_T_sp">get_Th_sp</a>


function OutputCell = ...
    get_T_sp_multicore(InputCell)

% Load some data from IAPWS-95
coeffs = readIAPWS95data();

% Check for the host mineral. If you want to change it, run set_fi_mineral
[~, pressureMinimum] = get_fi_mineral();

if length(InputCell) < 4; InputCell(4) = pressureMinimum; end;

Th_inf = InputCell(1);
V = InputCell(2)*1e-18;
T_sp_to_calc = InputCell(3);
T_sp_working = InputCell(4);

tolerance = get_tolerance;

% The fit options. Usually the fit converges after less than 10 iterations,
% if it doesn't there will be no minimum. The GradObj-entry tells the fit
% routine to take the Jacobian into account.
TolX = 1e-13;
TolFun = 1e-15;
options = optimset('TolX',TolX,'TolFun',TolFun,'GradObj','on',...
    'Hessian','user-supplied','Algorithm','trust-region-reflective',...
    'Display','off','MaxIter',50);

% Some constants
rhoc = 322;
Tc = 647.096;
rc = 1.1808741e-8;
b = -0.625;
mu = 1.256;

rhoOverallInitial = liqvap_density(Th_inf)*1000;

% Check whether there can be a bubble at all in the inclusion

% Apply the volume correction
[reftemp, alpha_V] = expansion_coeff(T_sp_working);
rho_overall_at_T = rhoOverallInitial*((1-(reftemp-Th_inf+273.15)*alpha_V)/(1-(reftemp-T_sp_working+273.15)*alpha_V));
dm = rho_overall_at_T/rhoc;

% Calculate the surface tension
tau = Tc/T_sp_working;
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
minvars_corrected(1) = (1 - rho_overall_at_T/1000/liqvap_density(T_sp_working));
minvars_corrected(2) = liqvap_density_vapor(T_sp_working)/rhoc*1000;

if minvars_corrected(1) > 0
    % Minimise the helmholtz energy
    minvars_corrected = fmincon(helmholtz_function, minvars_corrected, [],[],[],[],[0 0],[1 Inf],[],options);
    radius_out_corrected = minvars_corrected(1);
else
    radius_out_corrected = NaN;
end

if isnan(radius_out_corrected) || ~isreal(radius_out_corrected) || radius_out_corrected <= 1e-9
    % No minimum found; most probably there is no bubble possible
    % in this volume at this temperature;
    OutputCell(1) = NaN;
    OutputCell(2) = NaN;

    return;
end;
    
% There was a bubble possible at the density maximum, so proceed to
% look for the T_sp.

% First thing in the iteration will be to change the direction, so
% don't worry about the minus sign in the next line.
T_sp_working = InputCell(4);
dir = -T_sp_to_calc;
if abs(T_sp_working - pressureMinimum) < 0.5
    step = 1.25*5;
else
    step = 1.25;
end
while (T_sp_to_calc == 1 && T_sp_working > Th_inf);
    % This happens if the users initial guess is higher than Th_inf.
    % Usually this occurs when another function calls this routine
    T_sp_working = T_sp_working + dir*step;
end;

previous_radius = radius_out_corrected;

while step > tolerance
    % We stop as soon as we reached the tolerance in temperature.

    % We crossed T_sp; make the step smaller and change direction.
    step = step/5;
    dir = -dir;

    while isnan(radius_out_corrected) == isnan(previous_radius)
        % As long as we don't cross T_sp, the above condition
        % will hold. We keep on heating (or cooling for retrograde)
        % the inclusion and check if there is still a bubble
        
        if (sign(T_sp_working + dir*step - pressureMinimum) ~= sign(T_sp_working - pressureMinimum)) ...
                && ((T_sp_working - pressureMinimum) ~= 0)
            % We will cross the density maximum. This is bad. So
            % change to the density maximum.
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
            while (T_sp_working + dir*step >= Th_inf && T_sp_to_calc == 1);
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
        
    end;
    % When we reach here, we crossed T_sp.

    % Save the radius for the next iteration.
    temp_save_radius = previous_radius;
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

OutputCell(1) = T_sp;
OutputCell(2) = r_sp;

% Done.
return;
