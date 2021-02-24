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
        
        function [objective, mom_con, supp_con] = peak_cons(obj,d)
            %PEAKCONS formulate support and measure constraints for peak
            %program at degree d
            %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;

            supp_con = [];       %support constraint     
            mass_init_sum = 0;   %mass of initial measure should be 1
            objective = 0;
            
           
            liou_con = cell(length(obj.locations), 1);
            obj_con = cell(length(obj.locations), 1);
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
            if isnumeric(mass_init_sum)
                mass_con = [];
            else
                % mass_init_sum == 0 eliminates the constraint when there
                % is one location. keep it in.
                mass_con = (mass_init_sum - 1 == 0);
            end
            
            loc_con = [];
            for i = 1:length(liou_con)
                %TODO: figure out the sign convention for dual variables
                %should the negative sign be there on liou_con?
                loc_con = [loc_con; -liou_con{i} == 0; obj_con{i}];
            end
                
            
                        
            mom_con = [loc_con; mass_con; zeno_con];

        end    
    
        function [sol, dual_rec] = peak_solve(obj, objective, mom_con,supp_con)
            %PEAK_SOLVE formulate and solve peak estimation program from
            %constraints and objective    

            mset('yalmip',true);
            mset(sdpsettings('solver', obj.solver));

            P = msdp(max(objective), mom_con, supp_con);

            sol = struct;
            [sol.status,sol.obj_rec, ~,sol.dual_rec]= msol(P);        

        end

        function obj = dual_process(obj, order, dual_rec)
            %DUAL_PROCESS dispatch the dual variables from solution to
            %locations and measures, turn the variables into nonnegative
            %functions along trajectories
            
            %v: coefficients of liouville
            %beta: coefficients of cost (if able)
            %alpha: dual of zeno gaps
            
            rec_eq = dual_rec{1};
            rec_ineq = dual_rec{2};
            
            gamma = rec_eq(end);
            
            
            liou_offset = 0;
            cost_con_offset = 0;
            for i = 1:length(obj.locations)
                
                loc_curr = obj.locations{i};
                
                %liouville
                nvar_curr = 1 + length(loc_curr.vars.x);
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
            
            for i = 1:length(obj.guards)
%                 obj.guards{i}.zeno_dual = rec_ineq(cost_con_offset + i);
                obj.guards{i}.dual_process(rec_ineq(cost_con_offset + i));
            end
            
            
        end
        
        function sol = peak(obj, order)
            %the main call, the full peak program at the target order
            
            d = 2*order;
            [objective, mom_con, supp_con] = obj.peak_cons(d);
            sol = obj.peak_solve(objective, mom_con,supp_con);
            obj.dual_process(order, sol.dual_rec);
            
        end
        
        %% Sampler
        
        function [supp_loc, supp_g] = supp_g_eval(obj, t, x, id)
            %find the support evaluation of locations and guards at index
            %(or source index) id
            
            g_mask = find(cellfun(@(g) g.src.id == id, obj.guards));
            N_guards = length(g_mask);
            
            supp_g = zeros(N_guards, 1);
            supp_loc(1) = obj.locations{id}.supp_eval(t, x);
            for j = 1:N_guards
                supp_g(j) = obj.guards{g_mask(j)}.supp_eval(t, x);
            end
        end
        
        function [event_eval, terminal, direction] = loc_event(obj, t, x, id)                   
            %event function for @ode15 or other solver
            Npt = size(x, 2);
            event_eval = zeros(Npt);
            
            
            for i = 1:Npt
                tcurr = t(:, i);               
                xcurr = x(:, i);               
                
                [supp_loc, supp_g] = supp_g_eval(obj, tcurr, xcurr, id);
                event_eval(i) = supp_loc && all(~supp_g);
                
            end
                        
            %stop integrating when the system falls outside support
            
            terminal = 1;
            direction = 0;                        
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
                [~, supp_g] = obj.supp_g_eval(t_curr, x_curr, id_curr);
                if any(supp_g)
                    %there is a jump to another guard
                    %xcurr is in the support of some guard
                    
                    %returns the first valid guard
                    [~, g_id_new] = max(supp_g); %maybe random?
                    
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
                    
                    if zeno_count(g_id_new) >= g_new.zeno_cap
                        %maximum number of jumps is met or exceeded
                        break
                    end
                else
                    %no available guards to jump to
                    %outside support, end of trajectory
                    break
                end
            end
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
            if isnumeric(init_sampler.init)
                %given sample points
                N = size(init_sampler.init, 2);
                out_sim_multi = cell(N, 1);
                
                for i = 1:N
                    
                    x0 = init_sampler.init(:, 2);
                    if length(init_sampler.loc) == 1
                        init_loc = init_sampler.loc;
                    else
                        init_loc = init_sampler.loc(i);
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
%             out_sim = cell(init_sampler.N);
        end
        
        
    end
    
end

