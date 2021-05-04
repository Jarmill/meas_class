classdef manager_interface
    %MANAGER_INTERFACE A generic class to pose the LMI in measures.
    %   applications include variations on peak estimation
    %   Detailed explanation goes here
    
    properties
        loc;        %locations      
        solver = 'mosek';
        sdp_settings = [];
    end
    
    methods
        function obj = manager_interface(loc)
            %MANAGER_INTERFACE Construct an instance of this class
            %   set the location
            obj.loc = loc;
            
            obj.sdp_settings = sdpsettings('solver', obj.solver, 'mosek.MSK_DPAR_BASIS_TOL_S', 1e-8, ...
                'mosek.MSK_DPAR_BASIS_TOL_X', 1e-8, 'mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED', 1e-9, ...
                'mosek.MSK_DPAR_INTPNT_TOL_PATH', 1e-6);
        end
        
        function [optimal, mom_out, corner] = recover(obj, tol)
            if nargin < 2
                tol = 5e-4;
            end
            [optimal, mom_out, corner] = obj.loc.recover(tol);
            
        end
         
        function sol = solve(obj, objective, mom_con,supp_con)
            %SOLVE formulate and solve peak estimation program from
            %constraints and objective    
            
            %TODO: incorporate minquery into maximin (minimax) formulation

            mset('yalmip',true);
            %make sure the solution is precise
            mset(obj.sdp_settings);
            % https://docs.mosek.com/9.2/pythonfusion/parameters.html
            
            

            P = msdp(max(objective), mom_con, supp_con);

            sol = struct;
            [sol.status,sol.obj_rec, ~,sol.dual_rec]= msol(P);        
        end  
        
        function sol = run(obj, order, Tmax)
            %RUN the main call, the full peak program at the target order
            
            if nargin < 3
                Tmax = 1;
            end
            
            d = 2*order;
            
            %formulate constraints
            [objective, mom_con, supp_con, len_dual] = obj.cons(d, Tmax);
                        
            %solve the program
            sol = obj.solve(objective, mom_con,supp_con);
            
            %process dual variables
            obj = obj.dual_process(d, sol.dual_rec, len_dual);
        end
    end
    
    
    
    methods(Abstract)
        cons(obj, d, Tmax);
        %formulate support and measure constraints for 
        %program at degree d
                       
        
        dual_process(obj, d, dual_rec, len_dual);
        %DUAL_PROCESS dispatch the dual variables from solution to
        %locations and measures, turn the variables into nonnegative
        %functions along trajectories
    end
    
end

