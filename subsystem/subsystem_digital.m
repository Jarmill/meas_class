classdef subsystem_digital < subsystem_interface
    %SUBSYSTEM_DIGITAL A subsystem x+=f(t, x, th, w) of a possibly 
    %uncertain discrete dynamical system in measure analysis    
    
    %TODO: Debug and determine if this class is correct
    
    methods
        %% Constructor
        function obj = subsystem_digital(loc_supp, f, sys_id, loc_id)
            %SUBSYSTEM_DIGITAL Construct a digital (possibly uncertain) 
            %subsystem, fill in information
            
            %process input
            if nargin < 3
                sys_id = 1;
            end
            
            if nargin < 4
                loc_id = [];
            end
            
            %superclass constructor
            obj@subsystem_interface(loc_supp, f, sys_id, loc_id);
                                     
        end
        
        
        %% measure definition
        function meas_new = meas_def(obj, suffix)           
            %declare a variable for each measure
            vars_new = struct('t', [], 'x', [], 'th', [], 'w', []);           
            varnames = fields(vars_new);
            for i = 1:length(varnames)
                curr_name = varnames{i};
                curr_var = obj.vars.(curr_name);
                
                if ~isempty(curr_var)
                    %declare a new variable
                    new_name = [curr_name, obj.prefix, suffix];
                    mpol(new_name, length(curr_var), 1);
                    %load the new variable into vars_new
                    vars_new.(curr_name) = eval(new_name);
                end
%                 obj.vars.(curr_var) = vars.(curr_var);
            end
            
           
                %create new support as well
                supp_ref = obj.supp.supp_sys_pack(obj.supp.X_sys);
                supp_new = subs_vars(supp_ref, obj.supp.get_vars(), ...
                                [vars_new.t; vars_new.x; vars_new.th; vars_new.w]);
                       
            %define the measure
            meas_new = meas_base(vars_new, supp_new);
        end
        
        %% Constraints        
        function Ay = cons_liou(obj, d)
            %CONS_LIOU Liouville Equation includes a comparision with 
            %forward operator by a pushforward (discrete systems only)
            
            Ay = obj.meas_occ.mom_push(d, obj.vars, obj.f) - ...
                              obj.meas_occ.mom_monom_proj(d);
        end
                      
        
        
        %% Dual Recovery  
        
        function obj = dual_process(obj, v)
            %DUAL_PROCESS store dual functions and compute nonnegative
            %functions for this digital subsystem
            pushforward = eval(v, obj.vars.x, obj.f);
            Lv = pushforward - v;
            obj.nn = -Lv;
        end
        
        %% Evaluation and Sampling
        %TODO: Implement this
    end
end

