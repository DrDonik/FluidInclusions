%get_T_bin 
% [T_bin, r_bin] = 
%                   get_T_bin(Th_inf, V, T_bin_to_calc)
% calculates the binodal temperature of an inclusion of a given size and a
% given nominal homogenization temperature. Th_inf is an array of
% nominal homogenization temperatures (in Kelvin), V is an array of
% volumes (in micrometers³) for which the binodal temperatures should be
% calculated, tolerance is the accuracy the binodal temperature should be
% calculated, and T_bin_to_calc is whether the user wants prograde (1) or
% retrograde (-1) temperatures. Only the first two arguments are necessary,
% tolerance and T_bin_to_calc are optional.
%
% The output will be the binodal temperatures in T_bin (in Kelvin) and the
% bubble radii at these temperatures in r_bin (in micrometers).

function [T_bin, r_bin] = ...
    get_T_bin(Th_inf, V, T_bin_to_calc)

% Check wether the user forgot to convert to Kelvin
if Th_inf(1) < 240 || Th_inf(end) < 240
    disp('I think you forgot to convert your temperatures to Kelvin. Aborting.');
    T_bin = NaN;
    r_bin = NaN;
    return;
end;

tolerance = get_tolerance;

% Change this value to 1 for prograde or to -1 for retrograde calculation
if nargin <= 2 || isempty(T_bin_to_calc); T_bin_to_calc = 1; end; %prograde

% Load some data from IAPWS-95
coeffs = readIAPWS95data();

% Check for the host mineral. If you want to change it, run set_fi_mineral
[~, pressureMinimum] = get_fi_mineral();

caller = dbstack; %Check whether the function was called from the commandline

% Initialise the output arrays to zero.
T_bin = zeros(length(V),length(Th_inf));
if nargout > 1; r_bin = zeros(length(V),length(Th_inf)); end;

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

% Walk through all the supplied volumes.
for V_ctr = 1:length(V);
    
    V_working = V(V_ctr);
    
    if numel(V) > 1
        % If there are more than one volumes to calculate, display which
        % one is calculated.
        disp(['Calculating for V = ', num2str(V_working*1e18), 'um^3'])
    end;
    
    % Walk through all the nominal homogenisation temperatures
    for Th_ctr = 1:length(Th_inf);
        
        Th_inf_working = Th_inf(Th_ctr);
        if numel(Th_inf) > 1 && ~isempty(caller)
            disp(['Calculating for Th_inf = ', num2str(Th_inf_working-273.15), '°C'])
        end;
                
        rhoOverallInitial = liqvap_density(Th_inf_working)*1000;

        % Check whether there can be a bubble at all in the inclusion
        T_bin_working = pressureMinimum;

        % Apply the volume correction
        [reftemp, alpha_V] = expansion_coeff(T_bin_working);
        rho_overall_at_T = rhoOverallInitial*((1-(reftemp-Th_inf_working+273.15)*alpha_V)/(1-(reftemp-T_bin_working+273.15)*alpha_V));
        dm = rho_overall_at_T/rhoc;
        
        % Calculate the surface tension
        tau = Tc/T_bin_working;
        stprime = rc/V_working^(1/3) * ((tau - 1)/tau)^mu * ((1 + b)*tau - b);

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
            %minvars_corrected = fminunc(helmholtz_function, minvars_corrected, options);
            minvars_corrected = fmincon(helmholtz_function, minvars_corrected, [],[],[],[],[0 0],[1 Inf],[],options);
            gm_out_corrected = minvars_corrected(1);
        else
            gm_out_corrected = NaN;
        end
        
        if isnan(gm_out_corrected) || ~isreal(gm_out_corrected) || gm_out_corrected <= 1e-9
            % No minimum found; most probably there is no bubble possible
            % in this volume at this temperature;
            
            continue;
        elseif helmholtz_function(minvars_corrected) > helmholtz_function([0, minvars_corrected(2)])
            % The minimum was found, but it is already in the metastable
            % region
            
            continue;
        end;

        % There was a bubble possible at the density maximum and it was in
        % the stable region; start to look for T_bin
        
        % First thing in the iteration will be to change the direction, so
        % don't worry about the minus sign in the next line.
        dir = -T_bin_to_calc;
        step = 1.25*5;

        previous_gm = gm_out_corrected;
                
        while step/5 >= tolerance
            % We stop as soon as we reache the tolerance in temperature.

            % We crossed T_bin; make the step smaller and change
            % direction.
            step = step/5;
            dir = -dir;
            
            while isnan(gm_out_corrected) == isnan(previous_gm)
                % As long as we don't cross T_bin, the above condition
                % will hold. We keep on heating (or cooling for retrograde)
                % the inclusion and check if there is still a stable bubble

                if (sign(T_bin_working + dir*step - pressureMinimum) ~= sign(T_bin_working - pressureMinimum)) && ((T_bin_working - pressureMinimum) ~= 0)
                    % We will cross the density maximum. This is bad. So
                    % change to the density maximum.
                    T_bin_working = pressureMinimum;
                else
                    T_bin_working = T_bin_working + dir*step;
                    while (T_bin_working >= Th_inf_working && T_bin_to_calc == 1);
                        % Make sure we don't cross Th_inf, since we know we
                        % will not find a bubble beyond that temperature.
                        T_bin_working = T_bin_working - dir*step;
                        step = step/5;
                        T_bin_working = T_bin_working + dir*step;
                    end;
                end;
                
                % Apply the volume correction
                [reftemp, alpha_V] = expansion_coeff(T_bin_working);
                rho_overall_at_T = rhoOverallInitial*((1-(reftemp-Th_inf_working+273.15)*alpha_V)/(1-(reftemp-T_bin_working+273.15)*alpha_V));
                dm = rho_overall_at_T/rhoc;
                
                % Calculate the surface tension
                tau = Tc/T_bin_working;
                stprime = rc/V_working^(1/3) * ((tau - 1)/tau)^mu * ((1 + b)*tau - b);
                
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
                    %minvars_corrected = fminunc(helmholtz_function, minvars_corrected, options);
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

            end;
            % When we reach here, we crossed T_bin.
            
            % Save the radius for the next iteration.
            if nargout == 2; temp_save_gm = previous_gm; end;
            previous_gm = gm_out_corrected;
            
        end;
        
        % We found a T_bin (or gave up). Right now we have to find out
        % whether we're on the NaN-side of T_bin or on the side with a
        % stable bubble, and then we'll save the value for Th_obs accordingly.
        if isnan(gm_out_corrected)
            T_bin(V_ctr,Th_ctr) = T_bin_working-dir*step;
            if nargout == 2; r_bin(V_ctr,Th_ctr) = (3 * V_working * temp_save_gm / (4 * pi))^(1/3)*1e6; end;
        else
            T_bin(V_ctr,Th_ctr) = T_bin_working;
            if nargout == 2; r_bin(V_ctr,Th_ctr) = (3 * V_working * gm_out_corrected / (4 * pi))^(1/3)*1e6; end;
        end;

        % Finished for this Th_inf. Proceed with the next.
    end;

    if isempty(caller)
        if nargout == 1
            save('T_bin','Th_inf','V','T_bin');
        elseif nargout == 2
            save('T_bin','Th_inf','V','T_bin','r_bin','TolX','TolFun');
        end
    end;
    
% Finished for this volume. Proceed with the next.
end;


% For all entries that are still 0 there was no bubble possible, change
% them to NaN. (This makes plotting easier)
mask_T_bin = T_bin == 0;
T_bin(mask_T_bin) = NaN;

% And save the final output.
if isempty(caller)
    if nargout == 1
        save('T_bin','Th_inf','V','T_bin');
    elseif nargout == 2
        r_bin(mask_T_bin) = NaN;
        save('T_bin','Th_inf','V','T_bin','r_bin','TolX','TolFun');
    end;
end;

% Done.
return;
