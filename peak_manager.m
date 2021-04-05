classdef peak_manager
    %PEAK_MANAGER standard peak estimation manager    
    %only a single location 'loc'
    %   may reverse this
    
    properties
        loc;        
        solver = 'mosek';
    end
    
    methods
        function obj = peak_manager(loc_supp, f, objective)
            %PEAK_MANAGER Construct an instance of this class
            
            if nargin == 2
                objective = 0;
            end

            obj.loc = location(1, loc_supp, f, objective);
            
            obj.solver = 'mosek';
        end
        
        %% Formulating and solving program
        
        function [objective, mom_con, supp_con, len_liou, len_abscont] = peak_cons(obj,d, Tmax)
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

            supp_con = obj.loc.supp_con();       %support constraint     
            mass_init_sum = 0;   %mass of initial measure should be 1
                        
            %get moment constraints
            liou_con = -(obj.loc.liou_con(d));
            [objective, obj_con] = obj.loc.objective_con(d);
            [abscont_con, len_abscont] = obj.loc.abscont_con(d);

            %finalize moment constraints
            
            %mass of initial measure sums to one
            mass_init_con = (obj.loc.mass_init() - 1 == 0);

            %time independent: mass of sum of occupation measures are less
            %than Tmax. Implement/fix this?
            %TODO: get this done
                
            len_liou = length(liou_con);
            
            mom_con = [liou_con==0; abscont_con==0; obj_con; mass_init_con];  


        end                    
    
        function sol = peak_solve(obj, objective, mom_con,supp_con)
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
                        
            
            %TODO: correct the offsets for abscont and liouville
            %TODO: only one step through the for loop
            
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
            
        end
        
        function sol = peak(obj, order, Tmax)
            %the main call, the full peak program at the target order
            
            if nargin < 3
                Tmax = 1;
            end
            
            d = 2*order;
            [objective, mom_con, supp_con, len_liou, len_abscont] = obj.peak_cons(d, Tmax);
            
            
            sol = obj.peak_solve(objective, mom_con,supp_con);
            
%             gamma_ind =  length(mom_con) - length(obj.guards);
            gamma = sol.dual_rec{1}(len_liou+1);
%             obj.dual_process(order, sol.dual_rec, gamma);
            
            %TODO: dual process the abscont zeta functions 
        end       
        
        %% Sampler       
        
        function [event_eval, terminal, direction] = loc_event(obj, t, x, id)                   
            %event function for @ode15 or other solver
            Npt = size(x, 2);
            event_eval = zeros(Npt);
            
            
            for i = 1:Npt
                tcurr = t(:, i);               
                xcurr = x(:, i);                               
                event_eval = obj.loc.supp_eval(t, x);
            end
                        
            %stop integrating when the system falls outside support
            
            terminal = 1;                  
            direction = -1;   % negative direction
        end 
       
        
        function out_sim = sample_traj(obj, t0, x0, id0, Tmax)
            %SAMPLE_TRAJ Sample a single trajectory starting at (t0, x0) in
            %a specific location. Track the trajectory as it moves through
            %locations
            
            %OUTPUT:
            %out_sim is a struct holding the simulation output: time,
            %state, objective, and nonnegative functions from the dual
            %solution of SDP.
            
            %TODO: a simpler wrapper for sample_traj_loc
            
            t_curr = t0;
            x_curr = x0;

            out_sim.sim = {};

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
        
        function [out_sim_multi] = sample_traj_multi(obj, init_sampler, Tmax)
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
            
            %TODO: include uncertainty terms
            
            if isnumeric(init_sampler.init)
                %given sample points
                N = size(init_sampler.init, 2);
                out_sim_multi = cell(N, 1);               
                %parallel code requires splitting off separate objects
                for i = 1:N                    
                    x0 = init_sampler.init(2:end, i);
                    out_sim_multi{i} = obj.sample_traj(0, x0, 1, Tmax);
                end
                
            else
                %random sample.
                N = init_sampler.N;
                out_sim_multi = cell(N, 1);
                for i = 1:N                    
                    x0 = init_sampler.x();                    
                    out_sim_multi{i} = obj.sample_traj(0, x0, 1, Tmax);
                end
            end
            
        end
        
        
        %% Plotter
        
        %TODO: modify plotter code for uncertainty
        %very likely write separate class
        
%         % Plotter                
%         function plot_nonneg_loc(obj,osd)
%             FS_title = 14;
%             FS_axis = 12;
%             
%             % Locations
%             setup
%             nonneg_title = {'Initial Value', 'Decrease in Value', 'Cost Proxy'};
%             i = 0;
%             figure(50+i)
%             clf
%             ax_loc = obj.nonneg_axis_str(i);
%             for k = 1:3
%                 subplot(3, 1, k)
%                 hold on
%                 xlabel('time', 'FontSize', FS_axis)
%                 title(['Loc ', num2str(i), ' ', nonneg_title{k}], 'FontSize', FS_title)
%                 ylabel(ax_loc{k}, 'interpreter', 'latex', 'FontSize', FS_axis);
%             end            
%             for j = 1:length(osd.locations{i})
%                 traj_curr = osd.locations{i}{j};
%                 for k = 1:3            
%                     subplot(3, 1, k)
%                     plot(traj_curr.t, traj_curr.nonneg(:, k), 'c')
%                 end
%             end   
%         end
% 
% 
%     function str_out = nonneg_axis_str(obj,loc)
%         locs = num2str(loc);
%         str_out{1} = ['$\gamma - v(x)$'];
%         str_out{2} = ['$-L_{f', locs, '} v_', locs, '(x)$'];
%         str_out{3} = ['$v_', locs, '(x) - p_', locs, '(x)$' ];
%     end
            
            
        
        
    end
    
end

