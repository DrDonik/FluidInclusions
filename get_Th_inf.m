%get_Th_inf  Calculate Th_inf from Th_obs and r_obs
%
% inclusionObject = 
%              get_Th_inf(Th_obs, r_obs, Th_obs_is_T_bin, T_obs, mineralNumber, Th_inf, V)
%
% Calculates the nominal homogenisation temperature and the volume of a
% fluid inclusion. Only Th_obs and r_obs are mandatory inputs.
%   Th_obs: Th_obs of the inclusion (in deg C)
%       NB: This can also be a retrograde Th_obs (but this will be slower)
%   r_obs: r_obs of the inclusion (in um)
%   Th_obs_is_T_bin: Whether Th_obs corresponds to T_bin (1 or 0)
%   T_obs: The temperature r_obs was measured at (in deg C)
%   mineralNumber: The mineralNumber of the inclusion
%   Th_inf: The initial guess for Th_inf (in deg C)
%       NB: If you provided a retrograde Th_obs, you will still need to
%       provide the PROgrade initial guess
%   V: The initial guess for V (in um^3)
%
% The inclusionObject will be a fully valid inclusion object matching your 
% desired values.
%

function inclusionObject = get_Th_inf(Th_obs, r_obs, Th_obs_is_T_bin, T_obs, mineralNumber, Th_inf, V)
    %% Preliminaries, set up all the necessary values, if they are not given

    if length(Th_obs) ~= length(r_obs); inclusionObject = []; return; end;

    if nargin < 6
        if nargin < 5
            [mineralNumber, T_pressureMinimum] = inclusion.set_fi_mineral();
            if nargin < 4
                if nargin < 3; Th_obs_is_T_bin = 0; end;
                T_obs = T_pressureMinimum - 273.15;
            end
        else
            [mineralNumber, T_pressureMinimum] = inclusion.set_fi_mineral(mineralNumber);
        end
        Th_inf = Th_obs + 2;
        V = ones(size(Th_obs))*1e7;
    elseif length(Th_inf) ~= length(Th_obs) || length(V) ~= length(Th_obs)
        inclusionObject = [];
        return
    else
        [mineralNumber, T_pressureMinimum] = inclusion.set_fi_mineral(mineralNumber);
    end
    
    if length(T_obs) ~= length(Th_obs);
        if length(T_obs) == 1;
            T_obs = repmat(T_obs,1,length(Th_obs));
        else
            inclusionObject = [];
            return
        end
    end

    T_pressureMinimum = T_pressureMinimum - 273.15;

    if Th_obs < T_pressureMinimum
        Th_obs_is_Th_inf_r = 1;
        if nargin < 6; Th_inf = 2*T_pressureMinimum - Th_obs + 2; end;
    else
        Th_obs_is_Th_inf_r = 0;
    end

    tolerance = 1e-4;
    
    %% Start of the main routine

    for Th_obs_ctr = length(Th_obs):-1:1

        root_pos = [r_obs(Th_obs_ctr); Th_obs(Th_obs_ctr)];

        iterationCounter = 0;

        inclusionObject(Th_obs_ctr) = inclusion(Th_inf(Th_obs_ctr), V(Th_obs_ctr), mineralNumber);
        [r_obs_calculated, Th_obs_calculated] = get_r_obs_and_Th_obs(inclusionObject(Th_obs_ctr), T_obs(Th_obs_ctr), Th_obs_is_T_bin, Th_obs_is_Th_inf_r);

        Th_inf_step = -1;
        V_step = -V(Th_obs_ctr)*1e-2;

        while iterationCounter < 12

            iterationCounter = iterationCounter + 1;

            Th_obs_calculated_old = Th_obs_calculated;
            r_obs_calculated_old = r_obs_calculated;

            % Sometimes it happens that we cross pressureMinimum. This is bad, 
            % so prevent it from happening
            while Th_inf(Th_obs_ctr) - Th_inf_step < T_pressureMinimum
                Th_inf_step = Th_inf_step/2;
            end

            % Sometimes it happens that the V_step is too big and the volume would
            % become negative. This is bad, so prevent it from happening
            while V(Th_obs_ctr) - V_step < 0
                V_step = V_step/2;
            end

            Th_obs_calculated = NaN; % This is for NaN-Catching
            r_obs_calculated = NaN;

            while isnan(Th_obs_calculated) || isnan(r_obs_calculated)

                % First, do the step in Th_inf-direction
                Th_inf(Th_obs_ctr) = Th_inf(Th_obs_ctr) - Th_inf_step;

                % This is for NaN-Catching. If the flower-boundary is crossed,
                % Th_obs and radius_out_corrected become NaN. If only one of them
                % becomes NaN, you requested strange things, such as a different
                % pressureMinimum.
                Th_obs_calculated = NaN;
                r_obs_calculated = NaN;
                while isnan(Th_obs_calculated) || isnan(r_obs_calculated)

                    inclusionObject(Th_obs_ctr) = inclusion(Th_inf(Th_obs_ctr), V(Th_obs_ctr), mineralNumber);
                    [r_obs_calculated, Th_obs_calculated] = get_r_obs_and_Th_obs(inclusionObject(Th_obs_ctr), T_obs(Th_obs_ctr), Th_obs_is_T_bin, Th_obs_is_Th_inf_r);
                    
                    if isnan(Th_obs_calculated)
                        % we probably crossed the flower boundary
                        Th_inf_step = Th_inf_step/2;
                        Th_inf(Th_obs_ctr) = Th_inf(Th_obs_ctr) + Th_inf_step;
                    end

                end

                % Calculate the derivative
                r_obs_grad_Th_inf = (r_obs_calculated_old - r_obs_calculated)/Th_inf_step;
                Th_obs_grad_Th_inf = (Th_obs_calculated_old - Th_obs_calculated)/Th_inf_step;
                % go back
                Th_inf(Th_obs_ctr) = Th_inf(Th_obs_ctr) + Th_inf_step;

                % Now do the step in V-direction
                V(Th_obs_ctr) = V(Th_obs_ctr) - V_step;

                % This is for NaN-Catching. See above for details
                Th_obs_calculated = NaN;
                r_obs_calculated = NaN;
                while isnan(Th_obs_calculated) || isnan(r_obs_calculated)

                    inclusionObject(Th_obs_ctr) = inclusion(Th_inf(Th_obs_ctr), V(Th_obs_ctr), mineralNumber);
                    [r_obs_calculated, Th_obs_calculated] = get_r_obs_and_Th_obs(inclusionObject(Th_obs_ctr), T_obs(Th_obs_ctr), Th_obs_is_T_bin, Th_obs_is_Th_inf_r);

                    if isnan(Th_obs_calculated);
                        % we probably crossed the flower boundary
                        V_step = V_step/2;
                        V(Th_obs_ctr) = V(Th_obs_ctr) + V_step;
                    end

                end

                % Calculate the derivative
                r_obs_grad_V = (r_obs_calculated_old - r_obs_calculated)/V_step;
                Th_obs_grad_V = (Th_obs_calculated_old - Th_obs_calculated)/V_step;

                % What we just did is calculate the Jacobian:
                Jc = [r_obs_grad_Th_inf r_obs_grad_V; Th_obs_grad_Th_inf Th_obs_grad_V];

                % Now go to the new position (We didn't change the volume back ...)
                Th_inf(Th_obs_ctr) = Th_inf(Th_obs_ctr) - Th_inf_step;

                inclusionObject(Th_obs_ctr) = inclusion(Th_inf(Th_obs_ctr), V(Th_obs_ctr), mineralNumber);
                [r_obs_calculated, Th_obs_calculated] = get_r_obs_and_Th_obs(inclusionObject(Th_obs_ctr), T_obs(Th_obs_ctr), Th_obs_is_T_bin, Th_obs_is_Th_inf_r);

                % If one of the two is NaN here, we have to redo the whole thing,
                % since the combination of the two new values doesn't work.
                if isnan(Th_obs_calculated)
                    % we probably crossed the flower boundary, so let's go back
                    % and try again
                    V(Th_obs_ctr) = V(Th_obs_ctr) + V_step;
                    V_step = V_step/2;
                    Th_inf(Th_obs_ctr) = Th_inf(Th_obs_ctr) + Th_inf_step;
                    Th_inf_step = Th_inf_step/2;
                end

            end

            % The vectorised version of the function looks as follows
            F_vec = [r_obs_calculated; Th_obs_calculated];

            % The next step should take us to a better guess and will be calculated
            % using Newton's method
            % The backslash here is a left multiplication of the inverse of the
            % Jacobian

            next_step_vec = Jc\(F_vec - root_pos);

            Th_inf_step = next_step_vec(1);
            V_step = next_step_vec(2);

            if sum(((F_vec - root_pos).*[10; 1]).^2) < tolerance
                disp(['Deviation of Th_obs and r_obs smaller than ', num2str(tolerance)]);
                break
            end

            if iterationCounter == 12
                error('Too many iterations');
            end

            if sum((next_step_vec.*[100; 1]).^2) < tolerance^2;
                disp(['Step size in Th_inf and V smaller than ', num2str(tolerance/100)]);
                break
            end

        end

        disp(['Th_obs = ', num2str(Th_obs(Th_obs_ctr)) ,'C, r_obs = ', num2str(r_obs(Th_obs_ctr)), 'um']);
        disp(['Th_inf = ', num2str(Th_inf(Th_obs_ctr)), 'C, V = ', num2str(V(Th_obs_ctr)), 'um^3']);
        disp('');

    end

    return
end

%% Helper function to return the two calculated values r_obs and Th_obs, 
% for a given combination of Th_inf and V

function [r_obs_calculated, Th_obs_calculated] = get_r_obs_and_Th_obs(inclusionObject, T_obs, Th_obs_is_T_bin, Th_obs_is_Th_inf_r)
    inclusionObject.T = T_obs;
    r_obs_calculated = inclusionObject.r(inclusionObject.T == T_obs);

    if Th_obs_is_T_bin
        if Th_obs_is_Th_inf_r
            Th_obs_calculated = inclusionObject.T_bin_r;
        else
            Th_obs_calculated = inclusionObject.T_bin;
        end
    else
        if Th_obs_is_Th_inf_r
            Th_obs_calculated = inclusionObject.T_sp_r;
        else
            Th_obs_calculated = inclusionObject.T_sp;
        end
    end

    return
end
