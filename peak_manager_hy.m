classdef peak_manager_hy < handle
    %PEAK_MANAGER_HY Hybrid peak estimation manager
    %   Detailed explanation goes here
    
    properties
        locations;
        guards;
                
        solver;                
    end
    
    methods
        function obj = peak_manager_hy(locations_in,guards_in)
            %PEAK_MANAGER_HY Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;

            if ~iscell(locations_in)
                obj.locations = {locations_in};
            else
                obj.locations = locations_in;
            end
            
            if ~iscell(guards_in)
                obj.guards = {guards_in};                        
            else
                obj.guards = guards_in;                        
            end           
            
            obj.solver = 'mosek';
        end
        
        %% Formulating and solving program
        
        function [objective, mom_con, supp_con, len_liou] = peak_cons(obj,d, Tmax)
            %PEAKCONS formulate support and measure constraints for peak
            %program at degree d
            %Input:
            %   d:      Monomials involved in relaxation (2*order)
            %   Tmax:   Maximum time (only when time-independent)
            %
            %Output:
            %   objective:  target to maximize  (@mom)
            %   mom_con:    moment constraints  (@momcon)
            %   supp_con:   support constraints (@supcon)
            %   len_liou:   number of liouville constraints (uint32)

            supp_con = [];       %support constraint     
            mass_init_sum = 0;   %mass of initial measure should be 1
            objective = 0;                                   
            
            liou_con = cell(length(obj.locations), 1);
            obj_con = cell(length(obj.locations), 1);
            
            mass_occ_sum = 0;
            %process the location measures
            for i = 1:length(obj.locations)
                loc_curr = obj.locations{i};
                
                %accumulate support constraints
                supp_con = [supp_con; loc_curr.supp_con()];

                %find moment constraints of current location
                [obj_curr, obj_con_curr] = loc_curr.objective_con(d);                                
                
                if ~isempty(obj_curr)
                    objective = objective + obj_curr;
                end
                obj_con{i} = obj_con_curr;
                
                liou_con{i} = loc_curr.liou_con(d);  
                                 
                %initial measure has mass 1
                if ~isempty(loc_curr.meas_init)
                    mass_init_sum = mass_init_sum + loc_curr.meas_init.mass();
                end
                
                if isempty(loc_curr.vars.t)
                    mass_occ_sum = mass_occ_sum + loc_curr.meas_occ.mass();
                end
            end
            
            %process the guards
            zeno_con = []; %zeno mass constraints
            
            for i = 1:length(obj.guards)
                g_curr = obj.guards{i};
                
                %support constraints
                supp_con = [supp_con; g_curr.supp];
                
                %zeno
                zeno_con = [zeno_con; g_curr.zeno_con()];
                
                
                %liouville constraints
                [mom_src, mom_dest] = g_curr.liou_reset(d);
                
                liou_con{g_curr.src.id}  = liou_con{g_curr.src.id}  + mom_src;
                liou_con{g_curr.dest.id} = liou_con{g_curr.dest.id} + mom_dest;
            end
            
            %finalize moment constraints
            
            %mass of initial measure sums to one
            if isnumeric(mass_init_sum)
                mass_init_con = [];
            else
                % mass_init_sum == 0 eliminates the constraint when there
                % is one location. keep it in.
                mass_init_con = (mass_init_sum - 1 == 0);
            end
            
            %TIME-INDEPENDENT mass of occupation measures are less than
            %Tmax. Only when t is not included as a variable;
            if isnumeric(mass_occ_sum) && (nargin > 2)
                
                %TODO: access Tmax
                mass_occ_con = [];
            else
                % mass_init_sum == 0 eliminates the constraint when there
                % is one location. keep it in.
                mass_occ_con = (mass_occ_sum - Tmax == 0);
            end
            
            
            loc_con = [];
            liou_con_all = [];
            obj_con_all = [];
            for i = 1:length(liou_con)
                %TODO: figure out the sign convention for dual variables
                %should the negative sign be there on liou_con?
%                 loc_con = [loc_con; -liou_con{i} == 0; obj_con{i}];
                liou_con_all = [liou_con_all; -liou_con{i} == 0];
                obj_con_all  = [obj_con_all; obj_con{i}];
            end
                
            len_liou = length(liou_con_all);
            
                        
            mom_con = [liou_con_all; mass_init_con; mass_occ_con; obj_con_all; zeno_con];

        end                    
    
        function [sol, dual_rec] = peak_solve(obj, objective, mom_con,supp_con)
            %PEAK_SOLVE formulate and solve peak estimation program from
            %constraints and objective    

            mset('yalmip',true);
            %make sure the solution is precise
            mset(sdpsettings('solver', obj.solver, 'mosek.MSK_DPAR_BASIS_TOL_S', 1e-8, ...
                'mosek.MSK_DPAR_BASIS_TOL_X', 1e-8, 'mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED', 1e-9, ...
                'mosek.MSK_DPAR_INTPNT_TOL_PATH', 1e-6));
            % https://docs.mosek.com/9.2/pythonfusion/parameters.html
            P = msdp(max(objective), mom_con, supp_con);

            sol = struct;
            [sol.status,sol.obj_rec, ~,sol.dual_rec]= msol(P);        
        end
        
        function s_out = mmat_corner(obj)
            %get the top corner of the moment matrix for all measures
%             s_out = struct('locations',  cell(length(obj.locations), 1), ...
%                            'guards', cell(length(obj.locations), 1));
                       
            s_out = {};
            s_out.locations = cell(length(obj.locations), 1);
            s_out.guards = cell(length(obj.guards), 1);
            
            for i = 1:length(obj.locations)
                s_out.locations{i} = obj.locations{i}.mmat_corner();
            end
            
            for i = 1:length(obj.guards)
                s_out.guards{i} = obj.guards{i}.mmat_corner();
            end
            
        end
        
        function obj = dual_process(obj, order, dual_rec, gamma)
            %DUAL_PROCESS dispatch the dual variables from solution to
            %locations and measures, turn the variables into nonnegative
            %functions along trajectories
            
            %v: coefficients of liouville
            %beta: coefficients of cost (if able)
            %alpha: dual of zeno gaps
            
            rec_eq = dual_rec{1};
            rec_ineq = dual_rec{2};
            
            if nargin < 3
                gamma = rec_eq(end);
            end
                        
            
            liou_offset = 0;
            cost_con_offset = 0;
            for i = 1:length(obj.locations)
                
                loc_curr = obj.locations{i};
                
                %liouville
                %time independent
                if isempty(loc_curr.vars.t)
                    nvar_curr = length(loc_curr.vars.x);
                else
                    nvar_curr = length(loc_curr.vars.x)+1;
                end
                liou_len_curr = nchoosek(nvar_curr + 2*order, 2*order);
                
                v_coeff = rec_eq((1:liou_len_curr) + liou_offset);
                                
                liou_offset = liou_offset + liou_len_curr;
                
                %maximin cost duals
                n_obj = length(loc_curr.objective);
                if n_obj > 1
                    beta_curr = rec_ineq((1:n_obj) + cost_con_offset);
                    cost_con_offset = cost_con_offset + n_obj;
                elseif isempty(loc_curr.objective)
                    beta_curr = [];
                else
                    beta_curr = 1;
                end
                
                loc_curr.dual_process(order, v_coeff, beta_curr, gamma);
                
            end
            
            %TODO: dual variable for Tmax
            
            for i = 1:length(obj.guards)
%                 obj.guards{i}.zeno_dual = rec_ineq(cost_con_offset + i);
                obj.guards{i}.dual_process(rec_ineq(cost_con_offset + i));
            end
            
            
        end
        
        function sol = peak(obj, order, Tmax)
            %the main call, the full peak program at the target order
            
            if nargin < 3
                Tmax = 0;
            end
            
            d = 2*order;
            [objective, mom_con, supp_con, len_liou] = obj.peak_cons(d, Tmax);
            
            
            sol = obj.peak_solve(objective, mom_con,supp_con);
            
%             gamma_ind =  length(mom_con) - length(obj.guards);
            gamma = sol.dual_rec{1}(len_liou+1);
            obj.dual_process(order, sol.dual_rec, gamma);
            
        end
        
        %% Recovery
        
        function [optimal, mom_out, corner] = recover(obj, tol)
            %RECOVER if top corner of the moment matrix is rank-1, then
            %return approximate optimizer
            
            if nargin < 2
                tol = 5e-4;
            end
            
            optimal = zeros(length(obj.locations), 1);
            mom_out = cell(length(obj.locations), 1);
            corner = cell(length(obj.locations), 1);
            for i = 1:length(obj.locations)
                [optimal(i), mom_out{i}, corner{i}] = obj.locations{i}.recover(tol);                 
            end            
        end
        
        %% Sampler
        
        function [supp_loc] = supp_loc_eval(obj, t, x)
            %find the support evaluation of locations and guards at index
            
            
            N_loc = length(obj.locations);
            supp_loc = zeros(N_loc, 1);
%             supp_loc(1) = obj.locations{id}.supp_eval(t, x);
            for j = 1:N_loc
                supp_loc(j) = obj.locations{j}.supp_eval(t, x);
            end
        end
        
        function [supp_loc, supp_g, possible_g] = supp_g_eval(obj, t, x, id)
            %find the support evaluation of locations and guards at index
            %(or source index) id
            
            g_mask = find(cellfun(@(g) g.src.id == id, obj.guards));
            N_guards = length(g_mask);
            
            supp_g = zeros(N_guards, 1);
            possible_g = [];
            supp_loc(1) = obj.locations{id}.supp_eval(t, x);
            for j = 1:N_guards
                supp_g(j) = obj.guards{g_mask(j)}.supp_eval(t, x);
                if supp_g(j)
                    possible_g = [possible_g; obj.guards{g_mask(j)}.id];
                end
            end
        end
        
        function [event_eval, terminal, direction] = loc_event(obj, t, x, id)                   
            %event function for @ode15 or other solver
            Npt = size(x, 2);
            event_eval = zeros(Npt);
            
            
            for i = 1:Npt
                tcurr = t(:, i);               
                xcurr = x(:, i);               
                
                %assume that the guards on are on the boundary of the
                %region
%                 [supp_loc, supp_g] = supp_g_eval(obj, tcurr, xcurr, id);
%                 event_eval(i) = supp_loc && all(~supp_g);
                event_eval = obj.locations{id}.supp_eval(t, x);
            end
                        
            %stop integrating when the system falls outside support
            
            terminal = 1;
%             direction = 0;                        
            direction = -1;   % negative direction
        end 
       
        
        function out_sim = sample_traj(obj, t0, x0, id0, Tmax)
            %SAMPLE_TRAJ Sample a single trajectory starting at (t0, x0) in
            %a specific location. Track the trajectory as it moves through
            %locations
            %
            %OUTPUT:
            %out_sim is a struct holding the simulation output: time,
            %state, objective, and nonnegative functions from the dual
            %solution of SDP.
            
            t_curr = t0;
            x_curr = x0;
            id_curr = id0;
%             out_sim = struct('sim', {}, 'jump', {});
            out_sim.sim = {};
            out_sim.jump = {};
            zeno_count = zeros(length(obj.guards), 1);
            while (t_curr < Tmax)
                %sample within location
                loc_curr = obj.locations{id_curr};
                event_curr = @(t, x) obj.loc_event(t, x, id_curr); 
                out_sim_curr = loc_curr.sample_traj_loc(t_curr, x_curr, Tmax, event_curr);
                
                %store trajectory in location
                out_sim.sim{end+1} = out_sim_curr;
                
                t_curr = out_sim_curr.t(end);
                x_curr = out_sim_curr.x(end, :)';
                
                %figure out jump
                [~, supp_g, poss_g] = obj.supp_g_eval(t_curr, x_curr, id_curr);
                if any(supp_g)
                    %there is a jump to another guard
                    %xcurr is in the support of some guard
                    
                    %returns the first valid guard
%                     [~, g_id_new] = max(supp_g); %maybe random?
                    g_id_new = poss_g(1);
                    
                    zeno_count(g_id_new) = zeno_count(g_id_new) + 1;
                    g_new = obj.guards{g_id_new};
                    Rx = g_new.reset_eval(x_curr);
                    
                    %track the nonnegativity  
                    jump_curr = struct('t', t_curr, 'x', x_curr, 'x_jump', ...
                        Rx, 'guard', g_id_new);
                    
                    if g_new.dual.solved
                        %TODO
                        jump_curr.nonneg = g_new.nonneg(t_curr, x_curr);
                    end
                    out_sim.jump{end+1} = jump_curr;
                    
                    
                    %complete the jump
                    id_curr = g_new.dest.id;
                    x_curr = Rx;
                    
                    if zeno_count(g_id_new) > g_new.zeno_cap
                        %maximum number of jumps is exceeded
                        break
                    end
                else
                    %no available guards to jump to
                    %outside support, end of trajectory
                    break
                end
            end
            out_sim.t_end = t_curr;
        end
        
        function [out_sim_multi, out_sim_deal] = sample_traj_multi(obj, init_sampler, Tmax)
            %SAMPLE_TRAJ_MULTI sample multiple trajectories through the 
            %sample_traj. 
            %
            %Input: 
            %init_sampler is a struct
            %   N:      number of points to sample
            %   init:   (id, x) initial location/state funciton
            %
            %OR
            %   x:      state (array)
            %   loc:    location (scalar or array)
            %
            %OUTPUT:
            %   out_sim_multi:  a cell indexed by trajectory sample with
            %                   fields corresponding to trajectory
            %                   locations and guard jumps
            %   out_sim_deal:   A struct with fields 'locations', and
            %                   'guards' containing all trajectories in
            %                   each domain
            
            if nargin < 2
                Tmax = 1;
                parallel = 0;
            elseif nargin < 3
                parallel = 0;
            end
            
            if isnumeric(init_sampler.init)
                %given sample points
                N = size(init_sampler.init, 2);
                out_sim_multi = cell(N, 1);               
                %parallel code requires splitting off separate objects
                for i = 1:N                    
                    x0 = init_sampler.init(2:end, i);
                    if length(init_sampler.loc) == 1
                        id0 = init_sampler.loc;
                    else
                        id0 = init_sampler.loc(i);
                    end
                    out_sim_multi{i} = obj.sample_traj(0, x0, id0, Tmax);
                end
                
            else
                %random sample.
                N = init_sampler.N;
                out_sim_multi = cell(N, 1);
                for i = 1:N                    
                    [id0, x0] = init_sampler.init();                    
                    out_sim_multi{i} = obj.sample_traj(0, x0, id0, Tmax);
                end
            end
            
            %now `deal' the samples into locations and fields            
            %matlab hackery to make cells of empty cells
            out_sim_deal = struct;
            out_sim_deal.locations = cellfun(@num2cell, cell(length(obj.locations), 1), 'UniformOutput', false);
            out_sim_deal.guards = cellfun(@num2cell, cell(length(obj.guards), 1), 'UniformOutput', false);
            
            
            
            for i = 1:N %every sampled trajectory
                for j = 1:length(out_sim_multi{i}.sim) 
                    %every location the trajectory visits
                    traj_curr = out_sim_multi{i}.sim{j};
                    loc_id = traj_curr.id;
                    out_sim_deal.locations{loc_id} = [out_sim_deal.locations{loc_id}; traj_curr];
                end                
                
                for j = 1:length(out_sim_multi{i}.jump)
                    jump_curr = out_sim_multi{i}.jump{j};
                    jump_id = jump_curr.guard;
                    out_sim_deal.guards{jump_id} = [out_sim_deal.guards{jump_id}; jump_curr];
                end
            end
        end
        
        
        %% Plotter        
        function plot_nonneg_jump(obj,osd)
            % PLOT_NONNEG plot the nonnegative functions along the
            % sampled trajectories
            %osd: out_sim_deal
            
            FS_title = 14;
            FS_axis = 12;
            
            Ng = length(obj.guards);
            Nl = length(obj.locations);
            
            figure(20)
            clf
            %plot the guards
            for g = 1:Ng
                subplot(Ng, 1, g)
                hold on
                xlabel('time', 'FontSize', FS_axis)
                id_src  = obj.guards{g}.src.id;
                id_dest = obj.guards{g}.dest.id;
                ylabel(['$v_', num2str(id_src), '(x) - v_', num2str(id_dest), ...
                        '(R_', num2str(g),'(x))$'],'interpreter', 'latex', 'FontSize', FS_axis)
                title(['Guard ', num2str(obj.guards{g}.id), ' Transition'], 'FontSize', FS_title)

                for j = 1:length(osd.guards{g})
                        j_curr = osd.guards{g}{j};            
                        stem(j_curr.t, j_curr.nonneg, 'c') 
                end
            end

        end
        
        function plot_nonneg_loc(obj,osd)
                    FS_title = 14;
            FS_axis = 12;
            
            Ng = length(obj.guards);
            Nl = length(obj.locations);
            %% Locations
            %setup
            nonneg_title = {'Initial Value', 'Decrease in Value', 'Cost Proxy'};
            for i = 1:Nl
                figure(50+i)
                clf
                ax_loc = obj.nonneg_axis_str(i);
                for k = 1:3
                    subplot(3, 1, k)
                    hold on
                    xlabel('time', 'FontSize', FS_axis)
                    title(['Loc ', num2str(i), ' ', nonneg_title{k}], 'FontSize', FS_title)
                    ylabel(ax_loc{k}, 'interpreter', 'latex', 'FontSize', FS_axis);
                end            
                for j = 1:length(osd.locations{i})
                    traj_curr = osd.locations{i}{j};
                    for k = 1:3            
                        subplot(3, 1, k)
                        plot(traj_curr.t, traj_curr.nonneg(:, k), 'c')
                    end
                end   
            end
        end


    function str_out = nonneg_axis_str(obj,loc)
        locs = num2str(loc);
        str_out{1} = ['$\gamma - v_', locs, '(x)$'];
        str_out{2} = ['$-L_{f', locs, '} v_', locs, '(x)$'];
        str_out{3} = ['$v_', locs, '(x) - p_', locs, '(x)$' ];
    end
            
            
        
        
    end
    
end

