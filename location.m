classdef location < handle
    %LOCATION A location (space) of a hybrid system
    %   includes descriptions of the space as well as measures
    
    properties
        %location number
        id; 
        
        meas_init;
        meas_term;
        meas_occ;
                
        vars = struct('t', [], 'x', []);
        supp;
        
        f;
        objective;
        
        %still need to deal with dynamics f
        %and cost p (with maximin)
        
    end
    
    methods
        function obj = location(id, vars, supp, supp0, f, objective)
            %Location Construct an instance of this class
            %   Detailed explanation goes here
            
            %fill in properties
            obj.id = id;
            obj.vars = vars;
            obj.supp = supp;
            
            obj.f = f;
            obj.objective = objective;
            
            obj.meas_init = meas_base(obj.var_def('0', supp0));
            obj.meas_term = meas_base(obj.var_def('p', supp));
            obj.meas_occ  = meas_base(obj.var_def('occ', supp));
                                        
        end
        
        function [vars_new] = var_def(obj, suffix, supp_old)
            %VAR_DEF create new variables 't[suffix]_id',
            %'x[suffix]_id
            
            tname = ['t', suffix, '_', num2str(obj.id)];
            xname = ['x', suffix, '_', num2str(obj.id)];
            
            mpol(tname, 1, 1);
            mpol(xname, length(obj.vars.x), 1);
            
            t_new = eval(tname);
            x_new = eval(xname);
            
            vars_new= struct('t', t_new, 'x', x_new);
            
            if nargin == 3
                vars_new.supp = subs(supp_old, [obj.vars.t; obj.vars.x], ...
                                [t_new; x_new]);
            end
        end
        
        
        function cons = liou_con(obj, d, f)
            %LIOU_CON generate liouville constraint within location
            %
            %do not yet set this equal to zero (arithmetic operations not
            %are defined for constraints)
            if nargin == 2
                f = obj.f;
            end
            
            Ay_init =  obj.meas_init.mom_monom(d);
            Ay_term = -obj.meas_term.mom_monom(d);
            Ay_occ  =  obj.meas_occ.mom_lie(d, obj.vars, f);
            
            %TODO: Digital dynamics
            %Ay_occ  =  obj.meas_occ.mom_push(d, obj.vars, f);
            cons = Ay_init + Ay_term + Ay_occ;
        end
        
        function supp_con_out = supp_con(obj)
            supp_con_out = [obj.meas_init.supp;
                            obj.meas_term.supp;
                            obj.meas_occ.supp];
        end
        
        function [obj_max, obj_con] = objective_con(obj, d, objective)
            %OBJECTIVE_CON deal with the objective, which may be maximin
            if nargin == 2
                objective = obj.objective;
            end
            
            if length(objective) == 1            
                obj_max = mom(obj.meas_occ.var_sub(obj.vars, objective));
            
                obj_con = [];
            end
            
        end
%         
%         function [supp_con, mom_con, objective_out] = constraints(d)
%             %CONSTRAINTS generate support and measure constraints for this
%             %location 
%             
%             supp_con = obj.supp_con;
%             
%             mom_con = [];
%             
%             [objective_out, obj_con_curr] = obj.objective_con(d);
%             liou_con_curr = obj.liou_con(d);
%             
%             mom_con = [liou_con_curr; obj_con_curr];
%             
%         end
        
        function mass = mass_init(obj)
            mass = obj.meas_init.mass();
        end
        
        %something about processing dual_rec to get nonnegative functions
        
        
        %% overloads
        function e = isempty(obj)
            %is the support empty?
            %as in supp = []. The harder question would be 'does the basic
            %semialgebraic set formed by the constraints satisfy a
            %nullstellensatz?'
            e = isempty(obj.supp);
        end
    end
end

