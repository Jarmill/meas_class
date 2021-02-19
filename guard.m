classdef guard < meas_base
    %GUARD guard measure governing transitions
    %   Detailed explanation goes here
    
    properties
        id = [];
        src = [];
        dest = [];
        
        reset = [];
        
        zeno_cap = 100;
        
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
            
            obj.supp = subs(supp_old, [vars_old.t; vars_old.x], ...
                                [obj.vars.t; obj.vars.x]);
            obj.reset = subs(reset_old, [vars_old.t; vars_old.x], ...
                [obj.vars.t; obj.vars.x]);
            
            obj.meas = meas([obj.vars.t; obj.vars.x]);
        end
        
        function con = zeno_con(obj)
            %zeno execution constraint
            %at most (zeno_cap) transitions will occur on the guard
            con = [obj.mass() <= obj.zeno_cap];
        end
        
        function mom_out = reset_push(obj, d)
            v = obj.monom(d);
%             f_curr = obj.var_sub(vars_old, f_old);
            Rv = subs(v, [obj.vars.t; obj.vars.x], [obj.vars.t; obj.reset]);
            mom_out = mom(Rv);
        end
        
        function [mom_src, mom_dest] = liou_reset(obj, d)
            %LIOU_RESET Liouville expressions from the transition
            %mom_src:  from the source location
            %mom_dest: to the destination location
            mom_src = -obj.mom_monom(d);
            mom_dest = obj.reset_push(d);
            
        end
        
    end
end

