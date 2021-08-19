classdef sampler_hy
    %SAMPLER_HY Sample a hybrid system in multiple locations
    %   may include uncertainty
    
    properties
        guards;         %cell array of guard objects
        loc_smp;        %cell array of location sampler objects (e.g. sampler_uncertain)        
    end
    
    methods
        function obj = sampler_hy(sampler_in,guards_in)
            %SAMPLER_HY Construct an instance of this class
            %   Detailed explanation goes here
            obj.loc_smp = sampler_in;
            obj.guards = guards_in;
        end
        
        
      %% sampling trajectories
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
%                   event_curr = @(t, x) obj.loc_event(t, x, id_curr); 
%                 out_sim_curr = loc_curr.sample_traj_loc(t_curr, x_curr, Tmax, event_curr);
                out_sim_curr = obj.loc_smp{id_curr}.sample_traj(t_curr, x_curr, [], Tmax);
                
                loc_curr = obj.loc_smp{id_curr}.loc;
                
                if ~isempty(loc_curr.dual) && loc_curr.dual.solved
                    out_sim_curr.v = eval(loc_curr.dual.v, loc_curr.get_vars_end, [out_sim_curr.t/Tmax, out_sim_curr.x]');
                    out_sim_curr.nonneg = eval([loc_curr.dual.nn; loc_curr.sys{1}.dual.nn], loc_curr.get_vars, [out_sim_curr.t/Tmax, out_sim_curr.x, out_sim_curr.w]');
                end
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
                        jump_curr.nonneg = g_new.nonneg(t_curr/Tmax, x_curr);
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
        
   
%% support evaluation and events             
        function [supp_loc] = supp_loc_eval(obj, t, x)
            %find the support evaluation of locations and guards at index
            
            
            N_loc = length(obj.loc);
            supp_loc = zeros(N_loc, 1);
%             supp_loc(1) = obj.loc{id}.supp_eval(t, x);
            for j = 1:N_loc
                supp_loc(j) = obj.loc{j}.supp_eval(t, x);
            end
        end
        
        function [supp_loc, supp_g, possible_g] = supp_g_eval(obj, t, x, id)
            %find the support evaluation of locations and guards at index
            %(or source index) id
            
            g_mask = find(cellfun(@(g) g.src.id == id, obj.guards));
            N_guards = length(g_mask);
            
            supp_g = zeros(N_guards, 1);
            possible_g = [];
            supp_loc(1) = obj.loc_smp{id}.loc.supp_eval(t, x);
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
                event_eval = obj.loc{id}.supp_eval(t, x);
            end
                        
            %stop integrating when the system falls outside support
            
            terminal = 1;
%             direction = 0;                        
            direction = -1;   % negative direction
        end   
        
    function [out_sim_multi, out_sim_deal] = sample_traj_multi(obj, N, Tmax)
            %SAMPLE_TRAJ_MULTI sample multiple trajectories through the 
            %sample_traj. 
            %
            %Input: 
            %   N:      Number of points to sample
            %   Tmax:   Time horizon of trajectories
            %
            %OUTPUT:
            %   out_sim_multi:  a cell indexed by trajectory sample with
            %                   fields corresponding to trajectory
            %                   locations and guard jumps
            %   out_sim_deal:   A struct with fields 'locations', and
            %                   'guards' containing all trajectories in
            %                   each domain
            
            if nargin < 3
                Tmax = 1;
            end
            %TODO: no parallelization just yet. 
            
%             if nargin < 2
%                 Tmax = 1;
%                 parallel = 0;
%             elseif nargin < 3
%                 parallel = 0;
%             end
            
            %sample hybrid trajectories one at a time
            out_sim_multi = cell(N, 1); 
            init_loc_ind = cellfun(@(s) ~isempty(s.sampler.x), obj.loc_smp);
            elig_loc = (1:length(obj.loc_smp));
            elig_loc = elig_loc(init_loc_ind);
            for i = 1:N
                %uniformly select initial location out of list of eligible
                %locations. Set weights later.
%                 id0 = init_loc_ind(randi(length(init_loc_ind)));
                id0 = elig_loc(randi(length(elig_loc)));
                x0 = obj.loc_smp{id0}.sampler.x();
                
                out_sim_multi{i} = obj.sample_traj(0, x0, id0, Tmax);
    
            end

%             if isnumeric(init_sampler.init)
%                 %given sample points
%                 N = size(init_sampler.init, 2);
%                 out_sim_multi = cell(N, 1);               
%                 %parallel code requires splitting off separate objects
%                 for i = 1:N                    
%                     x0 = init_sampler.init(2:end, i);
%                     if length(init_sampler.loc) == 1
%                         id0 = init_sampler.loc;
%                     else
%                         id0 = init_sampler.loc(i);
%                     end
%                     out_sim_multi{i} = obj.sample_traj(0, x0, id0, Tmax);
%                 end
%                 
%             else
%                 %random sample.
%                 N = init_sampler.N;
%                 out_sim_multi = cell(N, 1);
%                 for i = 1:N                    
%                     [id0, x0] = init_sampler.init();                    
%                     out_sim_multi{i} = obj.sample_traj(0, x0, id0, Tmax);
%                 end
%             end
            
            %now `deal' the samples into locations and fields            
            %matlab hackery to make cells of empty cells
            out_sim_deal = struct;
            out_sim_deal.locations = cellfun(@num2cell, cell(length(obj.loc_smp), 1), 'UniformOutput', false);
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
      
    end
end

