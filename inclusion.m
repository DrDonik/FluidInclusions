%obj = inclusion(Th_inf, V, mineral_number)
%   The inclusion class is used to calculate, store and query all kinds
%   of information about inclusions. A new inclusion object is constructed
%   using inclusionObject = inclusion(Th_inf, V, mineral_number), the last
%   being an optional input (but asked during construction if omitted).
%   Th_inf has to be given in degree Celsius, V in um^3.
%   The new inclusionObject will contain nothing but the information
%   provided by the user at first, but is then be gradually filled with
%   information as soon as they are queried.
%
%   Properties:
%
%   Th_inf
%       The nominal homogenisation temperatures of the inclusionObject, in
%       degree C. This property is set during construction from user
%       input.
%   V
%       The volume of the inclusionObject, in um^3. This property is set
%       during construction from user input.
%   mineral
%       The mineral of the inclusionObject, i.e. the host type. This
%       property is set during construction from user input or query.
%   T_pressureMinimum
%       The temperature at which the internal pressure of the
%       inclusionObject is minimal, in degree C.
%   r_pressureMinimum
%       The radius of the vapour bubble at the pressure minimum, in um.
%   rho_overall
%       The mean density of the inclusionObject at Th_inf, in kg/m^3
%   T_sp, T_sp_r
%       The pro- and retrograde temperatures, respectively, at which 
%       the vapour bubble becomes unstable, in degree C.
%       This property is calculated as soon as it is queried.
%   r_sp, r_sp_r
%       The radiii of the vapour bubble at T_sp and T_sp_r, in um.
%       This property is calculated as soon as it is queried.
%   T_bin, T_bin_r
%       The pro- and retrograde temperatures, respectively, at which 
%       the vapour bubble becomes metastable, in degree C.
%       This property is calculated as soon as it is queried.
%   r_bin, r_bin_r
%       The radiii of the vapour bubble at T_bin and T_bin_r, in um.
%       This property is calculated as soon as it is queried.
%   Th_inf_r
%       The retrograde Th_inf, where the iso-Th-curve crosses the
%       saturation curve, in degree C.
%       This property is calculated as soon as it is queried.
%   flowerBoundary
%       The minimal Th_inf for this inclusion volume to still be able to
%       nucleate a bubble.
%   T
%       The "current" temperature(s) of the inclusionObject, in degree C.
%       This property is set using inclusion.T = values.
%   r
%       The bubble radii at temperature T, in um.
%       This property is calculated as soon as it is queried.
%   p_l, p_v, p_isoTh
%       The pressures of the liquid and vapour phase, respectively, at
%       the temperature(s) T, in Pascal. p_isoTh is the liquid pressure on
%       the iso-T_h curve (i.e. without bubble).
%       These properties are calculated as soon as it is queried.
%   rho_overall_at_T
%       The mean density of the inclusionObject at the temperature(s) T
%
%   Methods:
%
%   reset(inclusionObject)
%       This will clear all the calculated and hence stored values of
%       the limit temperatures and radii. Those will have to be
%       recalculated upon the next query. Useful if you changed
%       something in the underlying programs.
%   [inclusionObject, index, R] = errorTolerance(inclusionObject_array, Th_inf, Th_obs)
%       This will search through an inclusionObject array (consisting of
%       objects of the inclusion class) and return the entry that most
%       closely matches the Th_inf and Th_obs combination.
%   [Th_inf, V] = dimensions(inclusionObject_array)
%       This will return the span of Th_inf and V of an array of
%       inclusionObject.
%   inclusionObject_array = insert_Th_inf(inclusionObject_array, Th_inf)
%   inclusionObject_array = insert_V(inclusionObject_array, V)
%       The array of objects will be expanded so that it will contain
%       every possible combination of existing and inserted Th_inf and
%       V. That is, when adding a new, single Th_inf, the resulting 
%       object will contain a new full row with every combination of
%       that Th_inf with all preexisting V.
%   [inclusionObject, index, R] = find(inclusionObject_array, Th_inf, V)
%       This will search through an inclusionObject array (consisting of
%       objects of the inclusion class) and return the entry that most
%       closely matches the Th_inf and V combination.
%   [inclusionObject, index, R] = match(inclusionObject_array, Th_obs, r_obs, Th_obs_is_T_bin, T_obs)
%       This will search through an inclusionObject array (consisting of
%       objects of the inclusion class) and return the entry that most
%       closely matches the Th_obs and r_obs combination. If T_obs is not
%       given, r_obs will match r at the pressure minimum, otherwise it
%       will match r at T_obs.
%   inclusionObject = inclusion.get_Th_inf(Th_obs, r_obs, Th_obs_is_T_bin, T_obs, mineralNumber, Th_inf, V)
%       If you want an exact match (instead of an approximate one, see
%       above) use this function. It will calculate the inclusionObject
%       that fits the Th_obs and r_obs combination within an absolute
%       tolerance of 0.01 for both values. If T_obs is not given, r_obs
%       will match r at the pressure minimum, otherwise it will match r at
%       T_obs. Th_inf and V are start values for the iteration. You can use
%       match() (see above) to calculate them, it will speed up the process
%       by a few seconds.

classdef inclusion < hgsetget
    
    %% Properties
    
    properties (Dependent)
        Th_inf          % The nominal homogenization temperature of the inclusion, in degree C
        V               % The volume of the inclusion, in um^3
    end
    
    properties (SetAccess = immutable)  
        mineral         % The name of the mineral
    end
    
    properties (Dependent, SetAccess = protected)
        r_pressureMinimum   % The radius of the vapour bubble at the density maximum, in um
    end
    
    properties (SetAccess = immutable, Dependent, Hidden)
        T_pressureMinimum   % The temperature of the pressure minimum, in degree C
    end
    
    properties (SetAccess = immutable, Hidden)
        rho_overall     % The mean density of the liquid in the inclusion at Th_inf, in kg/m^3
    end
    
    properties (Dependent)
        T               % The "current" temperature(s) of the inclusion, in degree C
    end
    
    properties (Dependent, SetAccess = protected)
        r               % The radius of the vapour bubble at temperature(s) T, in um

        T_sp            % The temperature at which the bubble becomes unstable, in degree C
    end
	
    properties (Dependent, SetAccess = protected, Hidden)	
        r_sp            % The radius of the vapour bubble at T_sp, in um
        T_sp_r          % The retrograde temperature at which the bubble becomes unstable, in degree C
        r_sp_r          % The radius of the vapour bubble at T_sp_r, in um
        T_bin           % The temperature at which the bubble becomes metastable, in degree C
        r_bin           % The radius of the vapour bubble at T_bin, in um
        T_bin_r         % The retrograde temperature at which the bubble becomes metastable, in degree C
        r_bin_r         % The radius of the vapour bubble at T_bin_r, in um
        
        Th_inf_r        % The retrograde Th_inf, in degree C

        p_l             % The liquid pressure at temperature(s) T, in Pa
        p_v             % The vapour pressure at temperature(s) T, in Pa
        p_isoTh         % The liquid pressure without bubble at temperature(s) T, in Pa
        
        rho_overall_at_T        % The mean density of the inclusion at temperature(s) T, in kg/m^3
        flowerBoundary  % The minimum necessary Th_inf for this volume so that a bubble is possible
    end
    
    properties (SetAccess = immutable, GetAccess = protected)
        % The store_ properties are used to store once calculated values
        % of the different temperatures and radii.
        
        store_Th_inf    % This is to store Th_inf, in K
        store_V         % This is to store V, in m^3
        mineralNumber = 2;              % The host mineral
        store_T_pressureMinimum = 5.15+273.15;  % The temperature of maximum density, in K
    end
    
    properties (Access = protected)
        store_r_pressureMinimum      % The radius of the vapour bubble at the density maximum, in um
        store_T_sp      % The temperature at which the bubble becomes unstable, in K
        store_r_sp      % The radius of the vapour bubble at T_sp, in um
        store_T_sp_r    % The retrograde temperature at which the bubble becomes unstable, in K
        store_r_sp_r    % The radius of the vapour bubble at T_sp_r, in um
        store_T_bin     % The temperature at which the bubble becomes metastable, in K
        store_r_bin     % The radius of the vapour bubble at T_bin, in um
        store_T_bin_r   % The retrograde temperature at which the bubble becomes metastable, in K
        store_r_bin_r   % The radius of the vapour bubble at T_bin_r, in um

        store_Th_inf_r  % The retrograde Th_inf, in K
        
        store_T         % The "current" temperature(s) of the inclusion, in K
        store_r         % The radius of the vapour bubble at temperature(s) T, in um
        store_p_l       % The liquid pressure at temperature(s) T, in Pa
        store_p_v       % The vapour pressure at temperature(s) T, in Pa
        store_p_isoTh   % The liquid pressure without bubble at temperature(s) T, in Pa
        store_rho_overall_at_T       % The mean density of the inclusion at temperature(s) T, in kg/m^3
        store_flowerBoundary         % The minimum necessary Th_inf for this volume so that a bubble is possible
    end
    
    
    %% The constructor method
    methods
    
        function obj = inclusion(Th_inf, V, mineralNumber)

            if nargin == 0;
                Th_inf = 25+273.15;
                V = 1e6; 
                mineralNumber = 2;
            elseif nargin < 3
                mineralNumber = [];
            end
            
            if length(V) > 1 || length(Th_inf) > 1 || length(mineralNumber) > 1
                disp('One inclusion can only have one volume and homogenisation temperature');
                obj = [];
                return;
            end

            if isempty(mineralNumber)
                [obj.mineralNumber, obj.store_T_pressureMinimum, obj.mineral] = inclusion.set_fi_mineral();
            else
                [obj.mineralNumber, obj.store_T_pressureMinimum, obj.mineral] = inclusion.set_fi_mineral(mineralNumber);
            end
            
            obj.store_Th_inf = Th_inf + 273.15;
            obj.store_V = V/1e18;
            
            [~, obj.rho_overall] = saturationPressure(obj.store_Th_inf);
            
        end
        
    end
    
    %% Methods to display the values
    methods
        
        function pTplot(obj)
            % This will plot a pT-diagram with all the relevant curves
            
            % Make sure the inclusion contains all necessary T-values
            obj.T = [obj.Th_inf_r, obj.T_sp_r, obj.T_bin_r, obj.T_pressureMinimum, obj.T_bin, obj.T_sp, obj.Th_inf];
            
            figure
            plot(obj.T, obj.p_isoTh, 'k')
            hold on;
            plot(obj.T, obj.p_v, 'r')

            plot([obj.T_sp_r obj.T_sp_r], [obj.p_v(obj.T==obj.T_sp_r) obj.p_isoTh(obj.T==obj.T_sp_r)],'--k')

            plot(obj.T(find(obj.T == obj.T_sp_r,1):find(obj.T == obj.T_bin_r,1)-1), ...
                obj.p_l(find(obj.T == obj.T_sp_r,1):find(obj.T == obj.T_bin_r,1)-1),':b')
            
            plot([obj.T_bin_r obj.T_bin_r], [obj.p_v(obj.T==obj.T_bin_r) obj.p_isoTh(obj.T==obj.T_bin_r)],':k')

            plot(obj.T(find(obj.T == obj.T_bin_r,1):find(obj.T == obj.T_bin,1)-1), ...
                obj.p_l(find(obj.T == obj.T_bin_r,1):find(obj.T == obj.T_bin,1)-1),'b')
            
           plot([obj.T_bin obj.T_bin], [obj.p_v(obj.T==obj.T_bin) obj.p_isoTh(obj.T==obj.T_bin)],':k')

           plot(obj.T(find(obj.T == obj.T_bin,1):find(obj.T == obj.T_sp,1)-1), ...
                obj.p_l(find(obj.T == obj.T_bin,1):find(obj.T == obj.T_sp,1)-1),':b')
            
            plot([obj.T_sp obj.T_sp], [obj.p_v(obj.T==obj.T_sp) obj.p_isoTh(obj.T==obj.T_sp)],'--k')
            
            title(['Inclusion: ', obj.mineral, ', T_{h\infty} = ' num2str(obj.Th_inf), '\circC, V = ', num2str(obj.V), '\mum^3'])
            xlabel('T [degree C]')
            ylabel('p [Pa]')
        end
        
        function objProperties = getdisp(obj)
            % This will display all the properties that have been
            % calculated up to now. Use 'get(inclusionObject)' to use this.
            
            % set up the struct to display
            objProperties = struct;
            
            objProperties.Th_inf = obj.Th_inf;
            if ~isempty(obj.store_Th_inf_r); objProperties.Th_inf_r = obj.Th_inf_r; end;
            objProperties.V = obj.V;
            objProperties.mineral = obj.mineral;
            if ~isempty(obj.store_flowerBoundary); objProperties.flowerBoundary = obj.flowerBoundary; end;
            objProperties.T_pressureMinimum = obj.T_pressureMinimum;
            objProperties.r_pressureMinimum = obj.r_pressureMinimum;
            if ~isempty(obj.store_T_sp)
                objProperties.T_sp = obj.T_sp;
                objProperties.r_sp = obj.r_sp;
            end
            if ~isempty(obj.store_T_sp_r)
                objProperties.T_sp_r = obj.T_sp_r;
                objProperties.r_sp_r = obj.r_sp_r;
            end
            if ~isempty(obj.store_T_bin)
                objProperties.T_bin = obj.T_bin;
                objProperties.r_bin = obj.r_bin;
            end
            if ~isempty(obj.store_T_bin_r)
                objProperties.T_bin_r = obj.T_bin_r;
                objProperties.r_bin_r = obj.r_bin_r;
            end
            
            disp(objProperties);

        end
        
        
        function objProperties = getall(obj)
            % This will calculate all the properties and display them. Use
            % 'getall(inclusionObject)' to use this.
            
            % set up the struct to display
            objProperties = struct;
            
            objProperties.Th_inf = obj.Th_inf;
            objProperties.Th_inf_r = obj.Th_inf_r;
            objProperties.V = obj.V;
            objProperties.mineral = obj.mineral;
            objProperties.flowerBoundary = obj.flowerBoundary;
            objProperties.T_pressureMinimum = obj.T_pressureMinimum;
            objProperties.r_pressureMinimum = obj.r_pressureMinimum;
            objProperties.T_sp = obj.T_sp;
            objProperties.r_sp = obj.r_sp;
            objProperties.T_sp_r = obj.T_sp_r;
            objProperties.r_sp_r = obj.r_sp_r;
            objProperties.T_bin = obj.T_bin;
            objProperties.r_bin = obj.r_bin;
            objProperties.T_bin_r = obj.T_bin_r;
            objProperties.r_bin_r = obj.r_bin_r;
            
        end
        
    end
    
    %% Reset function
    methods
        
        function reset(obj)
            % This is used the clear all the calculated limit temperatures 
            % and radii. The values will then have to be calculated anew.
            % To reset T and the associated radii, assign T an empty array.
            obj.store_T_sp = [];
            obj.store_r_sp = [];
            obj.store_T_sp_r = [];
            obj.store_r_sp_r = [];
            obj.store_T_bin = [];
            obj.store_r_bin = [];
            obj.store_T_bin_r = [];
            obj.store_r_bin_r = [];
        end
        
    end
    
    %% Methods for array handling
    methods
        
        function [obj, index, R] = find(obj_array, Th_inf, V)
            % This will return the indices of the inclusion that most
            % closely matches Th_inf and V.            
            for V_ctr = size(obj_array,2):-1:1
                obj_Th_inf(:, V_ctr) = [obj_array(:,V_ctr).Th_inf];
                obj_V(:, V_ctr) = [obj_array(:,V_ctr).V];
            end
             
            [R, index] = min((obj_Th_inf(:)-Th_inf(:)).^2 + (obj_V(:)-V(:)).^2);
            
            obj = obj_array(index);

        end
        
        function [obj, index, R] = errorTolerance(obj_array, Th_inf, Th_obs, Th_obs_is_T_bin)
            % This will return the indices of the inclusion that most
            % closely matches Th_inf and Th_obs.        
            if nargin < 4
                Th_obs_is_T_bin = 0;
            end
            
            for V_ctr = size(obj_array,2):-1:1
                obj_Th_inf(:, V_ctr) = [obj_array(:,V_ctr).Th_inf];
                if Th_obs_is_T_bin
                    obj_Th_obs(:, V_ctr) = [obj_array(:,V_ctr).T_bin];
                else
                    obj_Th_obs(:, V_ctr) = [obj_array(:,V_ctr).T_sp];
                end
            end
             
            [R, index] = min((obj_Th_inf(:)-Th_inf(:)).^2 + (obj_Th_obs(:)-Th_obs(:)).^2);
            
            obj = obj_array(index);

        end
        
        function [obj, index, R] = match(obj_array, Th_obs, r_obs, Th_obs_is_T_bin, T_obs)
            % This will search through an array of inclusion objects and
            % return the entry that matches the requested combination of
            % Th_obs and r_obs most closely. Th_obs_is_T_bin decides
            % which temperature Th_obs represents and is 0 by default.
            % If T_obs is set r_obs will be compared to the radius at this
            % temperature instead of r_pressureMinimum.
            if nargin < 4
                Th_obs_is_T_bin = 0;
            elseif Th_obs_is_T_bin > 1;
                disp('Are you sure you want Th_obs_is_T_bin to be greater than one?')
                disp('I guess you wanted to tell me to look for a different T_obs.')
                disp('I will treat your query accordingly.')
                disp('(The proper way to call this function is ''match(obj_array, Th_obs, r_obs, Th_obs_is_T_bin, T_obs)'')')
                T_obs = Th_obs_is_T_bin;
                Th_obs_is_T_bin = 0;
            end

            if ~exist('T_obs','var')

                for V_ctr = size(obj_array,2):-1:1
                    r(:, V_ctr) = [obj_array(:,V_ctr).r_pressureMinimum];
                    if Th_obs_is_T_bin
                        obj_Th_obs(:, V_ctr) = [obj_array(:,V_ctr).T_bin];
                    else
                        obj_Th_obs(:, V_ctr) = [obj_array(:,V_ctr).T_sp];
                    end
                end
                
                for Th_obs_ctr = length(Th_obs):-1:1
                    [R(Th_obs_ctr), index(Th_obs_ctr)] = min((obj_Th_obs(:)-Th_obs(Th_obs_ctr)).^2 + (r(:)-r_obs(Th_obs_ctr)).^2);
                end
                
            else
                % Search for the entries that contain the requested temperature T_obs
 
                for V_ctr = size(obj_array,2):-1:1
                    for Th_inf_ctr = size(obj_array,1):-1:1
                        logical_indices(Th_inf_ctr, V_ctr) = ~isempty(find(abs(obj_array(Th_inf_ctr, V_ctr).T - T_obs) < 1e-2, 1));
                    end
                end

                if ~any(logical_indices)
                    disp('No entry matches the requested temperature')
                    decision = input('Do you want me to calculate all missing values? [y/N]: ', 's');
                    if strcmp(decision,'y') || strcmp(decision,'Y') || strcmp(decision,'yes') || strcmp(decision,'Yes') || strcmp(decision,'YES')
                        set(obj_array,'T',T_obs);
                        logical_indices = true(size(obj_array));
                    else
                        obj = [];
                        R = [];
                        index = [];
                        return;
                    end;
                elseif ~all(logical_indices)
                    disp('Not all entries in your array contain the requested temperature.')
                    decision = input('Do you want me to calculate all missing values? [y/N]: ', 's');
                    
                    if strcmp(decision,'y') || strcmp(decision,'Y') || strcmp(decision,'yes') || strcmp(decision,'Yes') || strcmp(decision,'YES')
                        set(obj_array,'T',T_obs);
                        logical_indices = true(size(obj_array));
                    else
                        disp('Your match will probably not be correct.')
                    end
                    
                end;
                
                sub_obj_array = obj_array(logical_indices);

                for r_ctr = length(sub_obj_array):-1:1
                    r(r_ctr) = sub_obj_array(r_ctr).r(find(abs(sub_obj_array(r_ctr).T - T_obs) < 1e-2, 1));
                    if Th_obs_is_T_bin
                        obj_Th_obs(r_ctr) = sub_obj_array(r_ctr).T_bin;
                    else
                        obj_Th_obs(r_ctr) = sub_obj_array(r_ctr).T_sp;                            
                    end
                end

                for Th_obs_ctr = length(Th_obs):-1:1
                    [R(Th_obs_ctr), sub_index(Th_obs_ctr)] = min((obj_Th_obs(:)-Th_obs(Th_obs_ctr)).^2 + (r(:)-r_obs(Th_obs_ctr)).^2);
                end
                                
                indices = find(logical_indices);
                index = indices(sub_index);
                
            end
            
            obj = obj_array(index);

        end
        
        function [Th_inf, V] = dimensions(obj_array)
            Th_inf = [obj_array(:,1).Th_inf];
            V = [obj_array(1,:).V];
            
            if length(unique(Th_inf)) ~= length(Th_inf) || length(unique(V)) ~= length(V)
                disp('Your array is formed improperly.')
                disp('Make sure it contains Th_inf in the first and V in the second dimension.')
                Th_inf = [];
                V = [];
            end;
            
        end
        
        function obj_array = insert_Th_inf(obj_array, Th_inf)
            % figure out the dimensions of the input array
            [Th_inf_old_array, V_array] = dimensions(obj_array);

            if isempty(Th_inf_old_array) || isempty(V_array)
                return
            end
            
            % sort all the (old and new) Th_inf values and make sure each
            % one is only once in the resulting array
            Th_inf_new_array = unique([Th_inf Th_inf_old_array]);
            
            if length(Th_inf_new_array) == length(Th_inf_old_array)
                return
            end
            
            % Look where the values were in the old array
            for Th_inf_ctr = length(Th_inf_new_array):-1:1
                old_Th_inf_index = find(Th_inf_old_array == Th_inf_new_array(Th_inf_ctr));
                if old_Th_inf_index
                    % copy the old row
                    obj_array(Th_inf_ctr, :) = obj_array(old_Th_inf_index, :);
                else
                    % create the new row
                    for V_ctr = length(V_array):-1:1
                        obj_array(Th_inf_ctr, V_ctr) = ...
                            inclusion(Th_inf_new_array(Th_inf_ctr), V_array(V_ctr), obj_array(1).mineralNumber);
                    end
                end
            end
        end
        
        function obj_array = insert_V(obj_array, V)
            % figure out the dimensions of the input array
            [Th_inf_array, V_old_array] = dimensions(obj_array);

            if isempty(Th_inf_array) || isempty(V_old_array)
                return
            end
            
            % sort all the (old and new) Th_inf values and make sure each
            % one is only once in the resulting array
            V_new_array = unique([V V_old_array]);
            
            if length(V_new_array) == length(V_old_array)
                return
            end
            
            % Look where the values were in the old array
            for V_ctr = length(V_new_array):-1:1
                old_V_index = find(V_old_array == V_new_array(V_ctr));
                if old_V_index
                    % copy the old row
                    obj_array(:, V_ctr) = obj_array(:, old_V_index);
                else
                    % create the new row
                    for Th_inf_ctr = length(Th_inf_array):-1:1
                        obj_array(Th_inf_ctr, V_ctr) = ...
                            inclusion(Th_inf_array(Th_inf_ctr), V_new_array(V_ctr), obj_array(1).mineralNumber);
                    end
                end
            end
        end
			
    end

    
    %% Set and Get methods for properties that will not change with temperature T
    methods
        
        function value = get.Th_inf(obj)
            value = obj.store_Th_inf - 273.15;
        end
        
        function value = get.Th_inf_r(obj)
            if isempty(obj.store_Th_inf_r)
                
                rho_diff = @(Th_inf_r_working) calc_rho_diff(obj, Th_inf_r_working);
                % The next line would be OK, if saturationDensity would
                % work below -36 degrees, which it doesn't at the moment.
                %Th_inf_r_working = [2*obj.store_T_pressureMinimum - obj.store_Th_inf-5 obj.store_T_pressureMinimum];
                Th_inf_r_working = [273.15-36 obj.store_T_pressureMinimum];
                
                obj.store_Th_inf_r = fzero(rho_diff, Th_inf_r_working);

            end
            value = obj.store_Th_inf_r - 273.15;
        end
        
        function value = get.V(obj)
            value = obj.store_V*1e18;
        end
        
        function value = get.T_pressureMinimum(obj)
            value = obj.store_T_pressureMinimum - 273.15;
        end
        
        function value = get.r_pressureMinimum(obj)
            if isempty(obj.store_r_pressureMinimum)
                obj.T = obj.T_pressureMinimum;
                obj.store_r_pressureMinimum = obj.r(obj.T == obj.T_pressureMinimum);
            end
            
            value = obj.store_r_pressureMinimum;
        end
        
        function value = get.flowerBoundary(obj)
            if isempty(obj.store_flowerBoundary)
                calculateFlowerBoundary(obj);
            end
            
            value = obj.store_flowerBoundary-273.15;
        end
        
        function value = get.T_sp(obj)
            if isempty(obj.store_T_sp)
                if isnan(obj.r_pressureMinimum)
                    obj.store_T_sp = NaN;
                    obj.store_r_sp = NaN;
                else
                    [obj.store_T_sp, obj.store_r_sp] = get_T_boundary(obj, 1, 1);
                end
            end
            
            value = obj.store_T_sp - 273.15;
        end

        function value = get.r_sp(obj)
            if isempty(obj.store_r_sp)
                if isnan(obj.r_pressureMinimum)
                    obj.store_T_sp = NaN;
                    obj.store_r_sp = NaN;
                else
                    [obj.store_T_sp, obj.store_r_sp] = get_T_boundary(obj, 1, 1);
                end
            end
            
            value = obj.store_r_sp;
        end

        function value = get.T_sp_r(obj)
            if isempty(obj.store_T_sp_r)
                if isnan(obj.r_pressureMinimum)
                    obj.store_T_sp_r = NaN;
                    obj.store_r_sp_r = NaN;
                else
                    [obj.store_T_sp_r, obj.store_r_sp_r] = get_T_boundary(obj, 1, 0);
                end
            end
            
            value = obj.store_T_sp_r - 273.15;
        end
        
        function value = get.r_sp_r(obj)
            if isempty(obj.store_r_sp_r)
                if isnan(obj.r_pressureMinimum)
                    obj.store_T_sp_r = NaN;
                    obj.store_r_sp_r = NaN;
                else
                    [obj.store_T_sp_r, obj.store_r_sp_r] = get_T_boundary(obj, 1, 0);
                end
            end
            
            value = obj.store_r_sp_r;
        end

        function value = get.T_bin(obj)
            if isempty(obj.store_T_bin)
                if isnan(obj.r_pressureMinimum)
                    obj.store_T_bin = NaN;
                    obj.store_r_bin = NaN;
                else
                    [obj.store_T_bin, obj.store_r_bin] = get_T_boundary(obj, 0, 1);
                end
            end

            value = obj.store_T_bin - 273.15;
        end

        function value = get.r_bin(obj)
            if isempty(obj.store_r_bin)
                if isnan(obj.r_pressureMinimum)
                    obj.store_T_bin = NaN;
                    obj.store_r_bin = NaN;
                else
                    [obj.store_T_bin, obj.store_r_bin] = get_T_boundary(obj, 0, 1);
                end
            end

            value = obj.store_r_bin;
        end

        function value = get.T_bin_r(obj)
            if isempty(obj.store_T_bin_r)
                if isnan(obj.r_pressureMinimum)
                    obj.store_T_bin_r = NaN;
                    obj.store_r_bin_r = NaN;
                else
                    [obj.store_T_bin_r, obj.store_r_bin_r] = get_T_boundary(obj, 0, 0);
                end
            end
            
            value = obj.store_T_bin_r - 273.15;
        end
    
        function value = get.r_bin_r(obj)
            if isempty(obj.store_r_bin_r)
                if isnan(obj.r_pressureMinimum)
                    obj.store_T_bin_r = NaN;
                    obj.store_r_bin_r = NaN;
                else
                    [obj.store_T_bin_r, obj.store_r_bin_r] = get_T_boundary(obj, 0, 0);
                end
            end
            
            value = obj.store_r_bin_r;
        end
        
    end


    %% Set and Get methods for properties that change with temperature T
    methods
        
        function set.T(obj, value)
            % Whenever the user assigns new "current" temperatures, we have
            % to reorder obj.store_r, obj.store_p_v, obj.store_p_l and 
            % obj.store_rho_overall_at_T since they will no longer match
            % with those temperatures. The old values of T will be kept, so
            % in fact you can just add values. To reset T, assign it an
            % empty array.
            
            % If an empty value is assigned, we reset T and all dependent values
            if isempty(value)
                obj.store_r = [];                                   
                obj.store_p_l = [];
                obj.store_p_v = [];
                obj.store_p_isoTh = [];
                obj.store_rho_overall_at_T = [];

                obj.store_T = [];
                return
            end;
            
            % We expect T to be input in C, convert it to K here
            value = value + 273.15;
            % We will need Th_inf_r, so trigger its calculation.
            %%obj.T_sp_r; obj.T_sp;
            
            temp_store_T = obj.store_T;
            
            value = unique([value temp_store_T]);
            
            obj.store_T = value;

            temp_store_r = obj.store_r;
            temp_store_p_l = obj.store_p_l;
            temp_store_p_v = obj.store_p_v;
            temp_store_p_isoTh = obj.store_p_isoTh;
            temp_store_rho_overall_at_T = obj.store_rho_overall_at_T;
            
            obj.store_r = zeros(size(value));                                   
            obj.store_p_l = zeros(size(value));
            obj.store_p_v = zeros(size(value));
            obj.store_p_isoTh = zeros(size(value));
            obj.store_rho_overall_at_T = zeros(size(value));
            
            % save all the values that can be kept
            for T_ctr = 1:length(value)
            
                index = find(temp_store_T==value(T_ctr), 1);
                if index
                    obj.store_r(T_ctr) = temp_store_r(index);                                   
                    obj.store_p_l(T_ctr) = temp_store_p_l(index);
                    obj.store_p_v(T_ctr) = temp_store_p_v(index);
                    obj.store_p_isoTh(T_ctr) = temp_store_p_isoTh(index);
                    obj.store_rho_overall_at_T(T_ctr) = temp_store_rho_overall_at_T(index);
                elseif 1==0%% value(T_ctr) >= obj.store_T_sp || value(T_ctr) <= obj.store_T_sp_r 
                    obj.store_r(T_ctr) = NaN;                                   
                    obj.store_p_l(T_ctr) = 0;
                    obj.store_p_v(T_ctr) = NaN;
                    obj.store_p_isoTh(T_ctr) = 0;
                    obj.store_rho_overall_at_T(T_ctr) = 0;
                end
                
            end
            
        end
        
        function value = get.T(obj)
            value = obj.store_T - 273.15;
        end

        function value = get.r(obj)
            if ~isempty(obj.store_T) && (isempty(obj.store_r) || ~all(obj.store_r))
                get_r(obj);
            end
            
            value = obj.store_r;
        end
               
        function value = get.rho_overall_at_T(obj)
            if ~isempty(obj.store_T) && (isempty(obj.store_rho_overall_at_T) || ~all(obj.store_rho_overall_at_T))
                [reftemp, alpha_V] = expansion_coeff(obj, obj.store_T(obj.store_rho_overall_at_T == 0));

                obj.store_rho_overall_at_T(obj.store_rho_overall_at_T == 0) = ...
                   obj.rho_overall * ((1-(reftemp-obj.store_Th_inf)*alpha_V) ./ ...
                   (1-(reftemp-obj.store_T(obj.store_rho_overall_at_T == 0)).*alpha_V));
            end
            
            value = obj.store_rho_overall_at_T;
        end
        
        function value = get.p_l(obj)
            if ~isempty(obj.store_T) && (isempty(obj.store_p_l) || ~all(obj.store_p_l))
                partPressure(obj);
            end
            
            value = obj.store_p_l;
        end

        function value = get.p_v(obj)
            if ~isempty(obj.store_T) && (isempty(obj.store_p_v) || ~all(obj.store_p_v))
                partPressure(obj);
           end
            
            value = obj.store_p_v;
        end
        
        function value = get.p_isoTh(obj)
            if ~isempty(obj.store_T) && (isempty(obj.store_p_isoTh) || ~all(obj.store_p_isoTh))
                partPressure(obj);
            end
            
            value = obj.store_p_isoTh;
        end

		
    end
    
    %% Helper methods
    methods (Access = protected)
        
        function rho_diff = calc_rho_diff(obj, Th_inf_r_working)
            
            [reftemp, alpha_V] = expansion_coeff(obj, Th_inf_r_working);

            rho_overall_at_Th_inf_r_working = ...
                obj.rho_overall * ((1-(reftemp-obj.store_Th_inf)*alpha_V) ./ ...
                (1-(reftemp-Th_inf_r_working).*alpha_V));

            [~, rho_sat] = saturationPressure(Th_inf_r_working);
           
            rho_diff = rho_sat - rho_overall_at_Th_inf_r_working;
            return
        end
        
    end
	
    %% Methods saved in files in the @inclusion folder
    methods (Access = protected)
        
        [reftemp, alpha_V] = expansion_coeff(obj, T)
        get_r(obj)
        calculateFlowerBoundary(obj)
        [T_boundary, r_boundary] = get_T_boundary(obj, calc_sp_boundary, calc_prograde_boundary)
        partPressure(obj)

    end
    
    %% static methods
    methods (Static)
		
        obj = get_Th_inf(Th_obs, r_obs, Th_obs_is_T_bin, T_obs, mineralNumber, Th_inf, V)
		
    end
	
    %% static helper methods
    methods (Static, Access = protected)
		
        coeffs = readIAPWS95data()
        [mineralNumber, pressureMinimum, mineral] = set_fi_mineral(mineralNumber)
        sigma = surface_tension(T)
                
    end
        
end
