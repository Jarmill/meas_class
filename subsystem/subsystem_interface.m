classdef subsystem_interface < handle
    %SUBSYSTEM_INTEFACE A generic subsystem x'=f(t, x,...)
    
    
    properties
        vars;
        
        %measures
        meas_occ = [];  %one occupation measure per subsystem
        
        f = [];         %dynamics
        
        %variables and identification
        loc_id = [];
        sys_id = [];
        prefix = '';                
        
        dual = struct('v', 0, 'Lv', 0, 'nn', 0);
        
        supp;       %support of measures
        supp_sys;   %support of only x
        
        meas_type = @meas_base;
    end
    
    properties(Access = protected)
        %use private properties for numeric function evaluations
        %object-oriented matlab is slow with public variables 
        %
        %TODO: check if protected is also slow
        f_;
        supp_sys_ = []; %support of state
        nn_ = 0;        %nonnegative function
    end
    
    methods
        %% Constructor
        function obj = subsystem_interface(loc_supp, f, sys_id, loc_id, meas_type)
            %SUBSYSTEM_INTERFACE Construct a subsystem, fill in information
            
            %process input
            if nargin < 3
                sys_id = 1;
            end
            
            if nargin < 4
                loc_id = [];
            end
            
            if nargin >= 5
                obj.meas_type = meas_type;
            end
            
            %identification and names
            obj.loc_id = loc_id;
            obj.sys_id = sys_id;
            obj.prefix = ['_', num2str(sys_id)];
            if ~isempty(loc_id)
                obj.prefix = ['_', num2str(loc_id), obj.prefix];
            end
            
            %support            

            obj.vars = loc_supp.vars;
            obj.supp_sys= loc_supp.get_X_sys_ind(sys_id);
            obj.supp_sys_ = obj.supp_sys;
            obj.supp = loc_supp.supp_sys_pack(obj.supp_sys);           
            %occupation measure definition
            obj.f = f;
            obj.f_ = f; 
            
            
            obj.meas_occ = obj.meas_def(obj.varnames, [obj.prefix,'_occ'], obj.supp);
%             obj.meas_occ  = obj.meas_def('occ');          
                                                                                
        end
        
        
        
        %% Getters
        
        function vars_out = get_vars(obj)
            %get variables in measure
            varnames = fields(obj.vars);
            vars_out = [];
            for i = 1:length(varnames)
                curr_var = varnames{i};
                vars_out = [vars_out; reshape(obj.vars.(curr_var), [], 1)];
            end
        end       
        
        function mass_out = mass_occ(obj)
            %mass of occupation measure (estimated total time)
            mass_out = obj.meas_occ.mass();
        end

        
        function abs_out = abscont_box(obj, d)
            %ABSCONT_BOX absolute continuity for box-inputs if present
            abs_out = [];
        end
              
        %% Getters
        function supp_all = get_supp(obj)
            %SUPP_ALL: get support set of all measures in subsystem
            
            supp_all = obj.meas_occ.supp;
        end        
        
        %% Measure Creation
         function meas_new = meas_def(obj, varnames, suffix, supp_ref)           
            %MEAS_DEF Define the measures in the collection
            %declare a variable for each measure (index ind in the union)

            vars_new = struct;
            old_stack =[];
            new_stack = [];

            for i = 1:length(varnames)
                curr_name = varnames{i};                
                
                if isfield(obj.vars, curr_name) 
                    if isempty(obj.vars.(curr_name))
                        vars_new.(curr_name) = [];
                    else
                        curr_var = obj.vars.(curr_name);
                        %declare a new variable
                        new_name = [curr_name, suffix];
                        mpol(new_name, size(curr_var, 1), size(curr_var, 2));
                        %load the new variable into vars_new
                        vars_new.(curr_name) = eval(new_name);

                        old_stack = [old_stack; reshape(obj.vars.(curr_name), [], 1)];
                        new_stack = [new_stack; reshape(vars_new.(curr_name), [], 1)];
                    end
                        
                end
            end
                       
            %get the support
            supp_new = subs_vars(supp_ref, old_stack, new_stack);
           
            
            %define the measure
            meas_new = obj.meas_type(vars_new, supp_new);
        end     
        
        %% Function Evaluation (for sampling)
        %TODO: Implement this
        %nonnegativity evaluation: include supports        
        
        %A function to evaluate dynamics f_
        function f_out = f_eval(obj, data)
            %data: [t, x, th, w, b] as required            
            f_out = eval(obj.f_, obj.get_vars(), data);
        end
        
        %A function to evaluate x in the support set supp_sys_
        function supp_out = supp_eval(obj, t, x)
            %is (t, x) in the support of the location?
                supp_out =  all(eval(obj.supp_sys_, obj.vars.x, x));
%             end
        end
        
        %An event function to detect when trajectories leave the set for
        %ODE solving
        function [event_eval, terminal, direction] = supp_event(obj, t, x)
            %event function for @ode15 or other solver
            Npt = size(x, 2);
            event_eval = zeros(1, Npt);
            for i = 1:Npt
                xcurr = x(:, i);
                tcurr = t(:, i);               

                event_eval(i) = obj.supp_eval(tcurr, xcurr);
            end
            
            %stop integrating when the system falls outside support
            
            terminal = 1;
            direction = 0;                        
        end
        
        %A function to evaluate nonnegative functions 
        function nn_out = nonneg_eval(obj, data)
            %data: [t, x, th, w] as required  
            nn_out = eval(obj.nn_, obj.get_vars(), data);
        end           
        
        
        
    end
    
    %% Overloads by inheritence
    methods(Abstract)
        
%         meas_new = meas_def(obj, suffix)           
        %declare a variable for each measure
        
        Ay = cons_liou(obj, d)
        %contribution of subsystem towards liouville constraints
        
        
        dual_process(obj, v, zeta)
        %DUAL_PROCESS store dual functions and compute nonnegative
        %functions for this subsystem
    end
end

