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
                obj.guards{i}.zeno_dual = rec_ineq(cost_con_offset + i);
            end
            
            
        end
        
        function sol = peak(obj, order)
            %the main call, the full peak program at the target order
            
            d = 2*order;
            [objective, mom_con, supp_con] = obj.peak_cons(d);
            sol = obj.peak_solve(objective, mom_con,supp_con);
            obj.dual_process(order, sol.dual_rec);
            
        end
        
    end
    
end

