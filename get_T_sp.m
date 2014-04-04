%get_T_sp
% [T_sp, r_sp] = 
%                 get_T_sp(Th_inf, V, T_sp_to_calc, initial_guess_T_sp, slowmo)
% calculates the spinodal temperature of an inclusion of a given size and a
% given nominal homogenization temperature. Th_inf is an array of
% nominal homogenization temperatures (in K), V is an array of
% volumes (in micrometers^3) for which the spinodal temperatures should be
% calculated, T_sp_to_calc is whether the user wants prograde (1) or
% retrograde (-1) temperatures, and initial_guess_T_sp gives an initial
% estimate for T_sp (in K). Only the first two arguments are necessary, 
% T_sp_to_calc and the initial guess are optional.
%
% This program will reformulate the statements and pass it on to 
% <a href="matlab: help get_T_sp_multicore">get_Th_sp_multicore</a>.
% You cann call the old, slow routine by setting the 5th input argument, slowmo, to 1.
%
% The output will be the spinodal temperatures in T_sp (in Kelvin) and the
% bubble radii at these spinodal temperatures in r_sp (in micrometers).
%
% See also <a href="matlab: help get_T_sp_multicore">get_T_sp_multicore</a>


function [T_sp, r_sp] = ...
    get_T_sp(Th_inf, V, T_sp_to_calc, initial_guess_T_sp, slowmo)

% Check wether the user forgot to convert to Kelvin
if Th_inf(1) < 240 || Th_inf(end) < 240
    disp('I think you forgot to convert your temperatures to Kelvin. Aborting.');
    T_sp = NaN;
    r_sp = NaN;
    return;
end;

tolerance = get_tolerance;

% Change this value to 1 for prograde or to -1 for retrograde calculation
if nargin <= 2 || isempty(T_sp_to_calc); T_sp_to_calc = ones(size(Th_inf)); end; %prograde

if nargin <= 5 || ~slowmo
    %disp (' ')
    %disp('You called the old version of get_T_sp.')
    %disp('I don''t think that''s what you wanted and will pass your arguments to get_T_sp_multicore.')
    if length(Th_inf) > 1
        disp('(Hint: you can start a few multicoreslaves now.)')
    end;
    
    if length(Th_inf) > length(V) && length(V) == 1
        outputCellArraySize = length(Th_inf);
        V = repmat(V,1,length(Th_inf)); 
        T_sp_to_calc = repmat(T_sp_to_calc,1,length(Th_inf));
    elseif length(V) > length(Th_inf) && length(Th_inf) == 1
        outputCellArraySize = length(V);
        Th_inf = repmat(Th_inf,1,length(V));
        T_sp_to_calc = repmat(T_sp_to_calc,1,length(V));
    elseif length(V) == length(Th_inf)
        outputCellArraySize = length(V);
    else
        disp('Sizes of V and Th_inf differ')
        return;
    end;
    
    inputCellArray = cell(1,outputCellArraySize);
    for ctr = 1:outputCellArraySize
        if nargin > 3
            inputCellArray{ctr} = [Th_inf(ctr), V(ctr), T_sp_to_calc(ctr), initial_guess_T_sp(ctr)];
        else
            inputCellArray{ctr} = [Th_inf(ctr), V(ctr), T_sp_to_calc(ctr)];
        end
    end
    outputCellArray = startmulticoremaster(@get_T_sp_multicore, inputCellArray);

    T_sp = zeros(size(Th_inf));
    r_sp = zeros(size(Th_inf));

    for ctr = 1:numel(outputCellArray)
        T_sp(ctr) = outputCellArray{ctr}(1);
        r_sp(ctr) = outputCellArray{ctr}(2);
    end

    return
end;

% Load some data from IAPWS-95
coeffs = readIAPWS95data();

% Check for the host mineral. If you want to change it, run set_fi_mineral
[~, pressureMinimum] = get_fi_mineral();

caller = dbstack; %Check whether the function was called from the commandline

% Initialise the output arrays to zero.
T_sp = zeros(length(V),length(Th_inf));
r_sp = zeros(length(V),length(Th_inf));

% The fit options. Usually the fit converges after less than 10 iterations,
% if it doesn't there will be no minimum. The GradObj-entry tells the fit
% routine to take the Jacobian into account.
TolX = 1e-13;
TolFun = 1e-15;
options = optimset('TolX',TolX,'TolFun',TolFun,'GradObj','on','Hessian','user-supplied','Display','off','MaxIter',50);

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
        T_sp_working = pressureMinimum;

        % Apply the volume correction
        [reftemp, alpha_V] = expansion_coeff(T_sp_working);
        rho_overall_at_T = rhoOverallInitial*((1-(reftemp-Th_inf_working+273.15)*alpha_V)/(1-(reftemp-T_sp_working+273.15)*alpha_V));
        dm = rho_overall_at_T/rhoc;
        
        % Calculate the surface tension
        tau = Tc/T_sp_working;
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
        minvars_corrected(1) = (1 - rho_overall_at_T/1000/liqvap_density(T_sp_working));
        minvars_corrected(2) = liqvap_density_vapor(T_sp_working)/rhoc*1000;
        
        if minvars_corrected(1) > 0
            % Minimise the helmholtz energy
            %minvars_corrected = fminunc(helmholtz_function, minvars_corrected, options);
            minvars_corrected = fmincon(helmholtz_function, minvars_corrected, [],[],[],[],[0 0],[1 Inf],[],options);
            radius_out_corrected = minvars_corrected(1);
        else
            radius_out_corrected = NaN;
        end
        
        if isnan(radius_out_corrected) || ~isreal(radius_out_corrected) || radius_out_corrected <= 1e-9
            % No minimum found; most probably there is no bubble possible
            % in this volume at this temperature;
            
            continue;
        end;

        % There was a bubble possible at the density maximum, so proceed to
        % look for the T_sp.
        
        % First thing in the iteration will be to change the direction, so
        % don't worry about the minus sign in the next line.
        dir = -T_sp_to_calc;
        step = 1.25*5;

        previous_radius = radius_out_corrected;
                
        while step/5 >= tolerance
            % We stop as soon as we reached the tolerance in temperature.

            % We crossed T_sp; make the step smaller and change
            % direction.
            step = step/5;
            dir = -dir;
            
            while isnan(radius_out_corrected) == isnan(previous_radius)
                % As long as we don't cross T_sp, the above condition
                % will hold. We keep on heating (or cooling for retrograde)
                % the inclusion and check if there is still a bubble

                if (sign(T_sp_working + dir*step - pressureMinimum) ~= sign(T_sp_working - pressureMinimum)) && ((T_sp_working - pressureMinimum) ~= 0)
                    % We will cross the density maximum. This is bad. So
                    % change to the density maximum.
                    T_sp_working = pressureMinimum;
                else
                    T_sp_working = T_sp_working + dir*step;
                    while (T_sp_working >= Th_inf_working && T_sp_to_calc == 1);
                        % Make sure we don't cross Th_inf, since we know we
                        % will not find a bubble beyond that temperature.
                        T_sp_working = T_sp_working - dir*step;
                        step = step/5;
                        T_sp_working = T_sp_working + dir*step;
                    end;
                end;
                
                % Apply the volume correction
                [reftemp, alpha_V] = expansion_coeff(T_sp_working);
                rho_overall_at_T = rhoOverallInitial*((1-(reftemp-Th_inf_working+273.15)*alpha_V)/(1-(reftemp-T_sp_working+273.15)*alpha_V));
                dm = rho_overall_at_T/rhoc;
                
                % Calculate the surface tension
                tau = Tc/T_sp_working;
                stprime = rc/V_working^(1/3) * (((tau - 1)/tau)^mu * ((1 + b)*tau - b));
                
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
                    %minvars_corrected = fminunc(helmholtz_function, minvars_corrected, options);
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
            if nargout == 2; temp_save_radius = previous_radius; end;
            previous_radius = radius_out_corrected;
            
        end;
        
        % We found a T_sp (or gave up). Right now we have to find out
        % whether we're on the NaN-side of T_sp or on the side with a
        % bubble, and then we'll save the value for T_sp accordingly.
        if isnan(radius_out_corrected)
            T_sp(V_ctr,Th_ctr) = T_sp_working-dir*step;
            if nargout == 2; r_sp(V_ctr,Th_ctr) = (3 * V_working * temp_save_radius / (4 * pi))^(1/3)*1e6; end;
        else
            T_sp(V_ctr,Th_ctr) = T_sp_working;
            if nargout == 2; r_sp(V_ctr,Th_ctr) = (3 * V_working * radius_out_corrected / (4 * pi))^(1/3)*1e6; end;
        end;

        % Finished for this Th_inf. Proceed with the next.
    end;

    if isempty(caller)
        if nargout == 1
            save('T_sp','Th_inf','V','T_sp');
        elseif nargout == 2
            save('T_sp','Th_inf','V','T_sp','r_sp','TolX','TolFun');
        end;
    end;

    
% Finished for this volume. Proceed with the next.
end;


% For all entries that are still 0 there was no bubble possible, change
% them to NaN. (This makes plotting easier)
mask_T_sp = T_sp == 0;
T_sp(mask_T_sp) = NaN;

% And save the final output.
if isempty(caller)
    if nargout == 1
        save('T_sp','Th_inf','V','T_sp');
    elseif nargout == 2
        r_sp(mask_T_sp) = NaN;
        save('T_sp','Th_inf','V','T_sp','r_sp','TolX','TolFun');
    end;
end;

% Done.
return;
