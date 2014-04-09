%get_Th_inf  Calculate Th_inf from Th_obs and r_obs
%
% inclusionObject = 
%              get_Th_inf(Th_obs, r_obs, Th_obs_is_T_bin, Th_inf, V)
%
% Calculates the nominal homogenisation temperature and the volume of a
% fluid inclusion
%   Th_obs: Th_obs of the inclusion (in C)
%   r_obs: r_obs of the inclusion (in um)
%   Th_obs_is_T_bin: Whether Th_obs corresponds to T_bin (1 or 0)
%   Th_inf: The initial guess for Th_inf (in C)
%   V: The initial guess for V (in um^3)
% All arguments are mandatory
%
% The inclusionObject will be a fully valid inclusion object matching your 
% desired values.
%

function inclusionObject = get_Th_inf(Th_obs, r_obs, Th_obs_is_T_bin, Th_inf, V)

if length(Th_obs) ~= length(r_obs); inclusionObject = []; return; end;

if nargin < 4
    if nargin < 3; Th_obs_is_T_bin = 0; end
    Th_inf = Th_obs + 5;
    V = ones(size(Th_obs))*1e7;
else
    if length(Th_inf) ~= length(Th_obs) || length(V) ~= length(Th_obs); inclusionObject = []; return; end;
end

[mineralNumber, pressureMinimum] = inclusion.set_fi_mineral();

tolerance = 1e-5;

for Th_obs_ctr = length(Th_obs):-1:1

    root_pos = [r_obs(Th_obs_ctr); Th_obs(Th_obs_ctr)];

    progress_counter = 0;

    inclusionObject(Th_obs_ctr) = inclusion(Th_inf(Th_obs_ctr), V(Th_obs_ctr), mineralNumber);

    if Th_obs_is_T_bin
        Th_obs_calculated = inclusionObject(Th_obs_ctr).T_bin;
    else
        Th_obs_calculated = inclusionObject(Th_obs_ctr).T_sp;
    end;
    radius_out_corrected = inclusionObject(Th_obs_ctr).r_pressureMinimum;

    Th_inf_step = -1;
    V_step = -V(Th_obs_ctr)/1e3;

    while abs(Th_obs_calculated - Th_obs) > tolerance || ...
            abs(radius_out_corrected - r_obs) > tolerance

        progress_counter = progress_counter+1;

        Th_obs_calculated_old = Th_obs_calculated;
        radius_out_corrected_old = radius_out_corrected;

        Th_obs_calculated = NaN; % This is for NaN-Catching
        radius_out_corrected = NaN;

        while isnan(Th_obs_calculated) || isnan(radius_out_corrected)

            % Sometimes it happens that we cross pressureMinimum. This is bad, so try
            % to recover:
            while Th_inf(Th_obs_ctr) - Th_inf_step < pressureMinimum
                Th_inf_step = Th_inf_step/2;
            end;

            % First, do the step in Th_inf-direction
            Th_inf(Th_obs_ctr) = Th_inf(Th_obs_ctr) - Th_inf_step;

            % This is for NaN-Catching. If the flower-boundary is crossed,
            % Th_obs and radius_out_corrected become NaN. If only one of them
            % becomes NaN, you requested strange things, such as a different
            % pressureMinimum.
            Th_obs_calculated = NaN;
            radius_out_corrected = NaN;
            while isnan(Th_obs_calculated) || isnan(radius_out_corrected)

                inclusionObject(Th_obs_ctr) = inclusion(Th_inf(Th_obs_ctr), V(Th_obs_ctr), mineralNumber);
                if Th_obs_is_T_bin
                    Th_obs_calculated = inclusionObject(Th_obs_ctr).T_bin;
                else
                    Th_obs_calculated = inclusionObject(Th_obs_ctr).T_sp;
                end;
                radius_out_corrected = inclusionObject(Th_obs_ctr).r_pressureMinimum;

                if isnan(Th_obs_calculated)
                    if sign(Th_inf_step) == 1;
                        % we probably crossed the flower boundary
                        Th_inf_step = Th_inf_step/2;
                        Th_inf(Th_obs_ctr) = Th_inf(Th_obs_ctr) + Th_inf_step;
                    else
                        % Something else went wrong.
                        keyboard;
                        %break;
                    end;
                elseif isnan(radius_out_corrected)
                    % You shouldn't end here, unless your PressureMinimum is not the
                    % original from set_fi_mineral
                    if Th_obs_calculated < pressureMinimum
                        % Told you. Let's try to recover.
                        Th_inf_step = Th_inf_step/2;
                        Th_inf(Th_obs_ctr) = Th_inf(Th_obs_ctr) + Th_inf_step;
                    else
                        % You definitely should never end here!
                        keyboard;
                        %break;
                    end;
                end;

            end;

            % Calculate the derivative
            r_grad_Th_inf = (radius_out_corrected_old - radius_out_corrected)/Th_inf_step;
            Th_obs_grad_Th_inf = (Th_obs_calculated_old - Th_obs_calculated)/Th_inf_step;
            % go back
            Th_inf(Th_obs_ctr) = Th_inf(Th_obs_ctr) + Th_inf_step;

            % Sometimes it happens that the V_step is too big and the volume would
            % become negative. This is bad, so prevent it from happening and define
            % an arbitrary step instead.
            while V(Th_obs_ctr) - V_step < 0
                V_step = V_step/2;
            end;

            % Now do the step in V-direction
            V(Th_obs_ctr) = V(Th_obs_ctr) - V_step;

            % This is for NaN-Catching. See above for details
            Th_obs_calculated = NaN;
            radius_out_corrected = NaN;
            while isnan(Th_obs_calculated) || isnan(radius_out_corrected)

                inclusionObject(Th_obs_ctr) = inclusion(Th_inf(Th_obs_ctr), V(Th_obs_ctr), mineralNumber);
                if Th_obs_is_T_bin
                    Th_obs_calculated = inclusionObject(Th_obs_ctr).T_bin;
                else
                    Th_obs_calculated = inclusionObject(Th_obs_ctr).T_sp;
                end;
                radius_out_corrected = inclusionObject(Th_obs_ctr).r_pressureMinimum;

                if isnan(Th_obs_calculated);
                    if sign(V_step) == 1;
                        % we probably crossed the flower boundary
                        V_step = V_step/2;
                        V(Th_obs_ctr) = V(Th_obs_ctr) + V_step;
                    else
                        % Something else went wrong.
                        keyboard;
                        %break;
                    end;
                elseif isnan(radius_out_corrected)
                    % You shouldn't end here, unless your PressureMinimum is not the
                    % original from set_fi_mineral
                    if Th_obs_calculated < pressureMinimum
                        % Told you. Let's try to recover.
                        V_step = V_step/2;
                        V(Th_obs_ctr) = V(Th_obs_ctr) + V_step;
                    else
                        % You definitely should never end here!
                        keyboard;
                        %break;
                    end;
                end;

            end;

            % Calculate the derivative
            r_grad_V = (radius_out_corrected_old - radius_out_corrected)/V_step;
            Th_obs_grad_V = (Th_obs_calculated_old - Th_obs_calculated)/V_step;

            % What we just did is calculate the Jacobian:
            Jc = [r_grad_Th_inf r_grad_V; Th_obs_grad_Th_inf Th_obs_grad_V];

            % Now go to the new position (We didn't change the colume back ...)
            Th_inf(Th_obs_ctr) = Th_inf(Th_obs_ctr) - Th_inf_step;

            inclusionObject(Th_obs_ctr) = inclusion(Th_inf(Th_obs_ctr), V(Th_obs_ctr), mineralNumber);
            if Th_obs_is_T_bin
                Th_obs_calculated = inclusionObject(Th_obs_ctr).T_bin;
            else
                Th_obs_calculated = inclusionObject(Th_obs_ctr).T_sp;
            end;
            radius_out_corrected = inclusionObject(Th_obs_ctr).r_pressureMinimum;

            % If one of the two is NaN here, we have to redo the whole thing,
            % since the combination of the two new values doesn't work.
            if isnan(Th_obs_calculated)
                % we probably crossed the flower boundary, so let's go back
                % and try again
                V(Th_obs_ctr) = V(Th_obs_ctr) + V_step;
                V_step = V_step/2;
                Th_inf(Th_obs_ctr) = Th_inf(Th_obs_ctr) + Th_inf_step;
                Th_inf_step = Th_inf_step/2;
            elseif isnan(radius_out_corrected);
                % You shouldn't end here, unless your PressureMinimum is not the
                % original from set_fi_mineral
                if Th_obs_calculated < pressureMinimum
                    % Told you. Let's try to recover.
                    Th_inf_step = Th_inf_step/2;
                    Th_inf(Th_obs_ctr) = Th_inf(Th_obs_ctr) + Th_inf_step;
                else
                    % You definitely should never end here! If you do, I have
                    % to give up.
                    keyboard;
                    %break;
                end;
            end;

        end;

        % The vectorised version of the function looks as follows
        F_vec = [radius_out_corrected; Th_obs_calculated];

        % The next step should take us to a better guess and will be calculated
        % using Newton's method
        % The backslash here is a left multiplication of the inverse of the
        % Jacobian

        next_step_vec = Jc\(F_vec - root_pos);

        if isnan(next_step_vec(1)) || isnan(next_step_vec(2));
            % We'll end here if JC doesn't have an inverse.
            % Maybe one of the steps is too small, then we don't have to make
            % that step. Let's figure it out:
            keyboard
            %Jc = [r_grad_Th_inf r_grad_V; Th_obs_grad_Th_inf Th_obs_grad_V];

            %next_step_vec = pinv(Jc)*(F_vec - root_pos);
        end;

        Th_inf_step = next_step_vec(1);
        V_step = next_step_vec(2);

        if abs(Th_inf_step/Th_inf(Th_obs_ctr)) < tolerance/100 && ...
                abs(V_step/V(Th_obs_ctr)) < tolerance/100;
            disp('Relative step size in Th_inf and relative step size in V smaller than 1% of tolerance');
            break;
        end

        if progress_counter >= 12
            disp('Too many iterations');
            break;
        end;

    end;
    
end;

return;