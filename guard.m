classdef guard < meas_base
    %GUARD guard measure governing transitions
    %   Detailed explanation goes here
    
    properties
        id = [];
        src = [];
        dest = [];
        
        reset = [];
        reset_identity = 0;
        
        zeno_cap = 2; %maximum number of transitions along guard
        dual = struct('zeno', 0, 'solved', 0);  %dual variable to zeno constraints
        
    end
    
    methods
        function obj = guard(id, vars_old, src,dest,supp_old, reset_old)
            %GUARD Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            
            obj@meas_base(vars_old, supp_old);
            obj.id = id;
            obj.src = src;
            obj.dest = dest;
%             obj.reset = reset_old;
            
            obj = obj.var_def('g', vars_old, supp_old, reset_old);
                        
        end
        
        function obj = set_zeno(zeno_new)
            obj.zeno_cap = zeno_new;
        end
        
        function obj = var_def(obj, suffix, vars_old, supp_old, reset_old)
            %VAR_DEF create new variables 't[suffix]_id',
            %'x[suffix]_id
            
            tname = ['t', suffix, '_', num2str(obj.id)];
            xname = ['x', suffix, '_', num2str(obj.id)];
            
            mpol(tname, 1, 1);
            mpol(xname, length(obj.vars.x), 1);
            
            obj.vars.t = eval(tname);
            obj.vars.x = eval(xname);
            
            obj.supp = subs_vars(supp_old, [vars_old.t; vars_old.x], ...
                                obj.get_vars());
            obj.reset = subs_vars(reset_old, [vars_old.t; vars_old.x], ...
                obj.get_vars());
            
            
            
            obj.meas = meas(obj.get_vars());
        end
        
        function con = zeno_con(obj)
            %zeno execution constraint
            %at most (zeno_cap) transitions will occur on the guard
            con = [obj.mass() <= obj.zeno_cap];
        end
        
        function mom_out = reset_push(obj, d)
            v = obj.monom(d);
%             f_curr = obj.var_sub(vars_old, f_old);
            Rv = subs(v, obj.get_vars(), [obj.vars.t; obj.reset]);
            mom_out = mom(Rv);
        end
        
        function [mom_src, mom_dest] = liou_reset(obj, d)
            %LIOU_RESET Liouville expressions from the transition
            %mom_src:  from the source location
            %mom_dest: to the destination location
            mom_src = -obj.mom_monom(d);
            mom_dest = obj.reset_push(d);
            
        end
        
        function reset_out = reset_eval(obj, x)
            %output of the reset map
            reset_out =  (eval(obj.reset, obj.vars.x, x));
        end
        
        function supp_out = supp_eval(obj, t, x)
            %is (t, x) in the support of the guard?
            supp_out =  all(eval(obj.supp, obj.get_vars(), [t; x]));
        end
        
        function nn_out  = nonneg(obj, t, x)
            %nonnegative dual function at this state transition
            %x comes from the space of src
            
            %check the sign convention
            %v should decrease along the jump
            
            Rx = eval(obj.reset, obj.vars.x, x);
            
            vsrc = obj.src.v_eval(t, x);
            vdest = obj.dest.v_eval(t, Rx);
            
            nn_out = vsrc - vdest - obj.dual.zeno;
            
        end
        
        function obj = dual_process(obj, zeno_dual)
            %store the dual variable
            obj.dual = struct('zeno', zeno_dual, 'solved', 1);
        end
                
        
    end
end

