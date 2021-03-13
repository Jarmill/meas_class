classdef location < handle
    %LOCATION A location (space) of a hybrid system
    %   includes descriptions of the space as well as measures
    
    properties
        %location number
        id; 
        
        %measures
        meas_init;
        meas_term;
        meas_occ;
        
        %variables
        vars = struct('t', [], 'x', []);
        cost_q = []; %multiple costs
        supp;   %support
        
        f;          %dynamics
        objective;
        
        dual = struct('v', 0, 'Lv', 0, 'beta', 1, 'gamma', 0, 'solved', 0);         
        %still need to deal with dynamics f
        %and cost p (with maximin)        
    end
    
    properties(Access = private)
        TIME_INDEP = 0;
    end
    
    methods
        function obj = location(id, vars, supp, supp0, f, objective)
            %Location Construct an instance of this class
            %   Detailed explanation goes here
            
            %fill in properties
            obj.id = id;
            obj.vars = vars;
            if ~isfield(vars, 't') || isempty(vars.t)
                obj.vars.t = []; %time-independent
                obj.TIME_INDEP = 1;            
            end
            obj.supp = supp;
            
            obj.f = f;            
            obj.objective = objective;
            
            if ~isempty(supp0)
                obj.meas_init = meas_base(obj.var_def('0', supp0));
            end
            
            if ~isempty(objective)
                obj.meas_term = meas_base(obj.var_def('p', supp));
            end
            obj.meas_occ  = meas_base(obj.var_def('occ', supp));                                            
        end
        
        function vars_out = get_vars(obj)
            %GET_VARS add more variables as necessary
            vars_out = [obj.vars.t; obj.vars.x];
        end
        
        function [vars_new] = var_def(obj, suffix, supp_old)
            %VAR_DEF create new variables 't[suffix]_id',
            %'x[suffix]_id
            
            if isempty(obj.vars.t)
                t_new = [];
            else
                tname = ['t', suffix, '_', num2str(obj.id)];                       
                mpol(tname, 1, 1);
                t_new = eval(tname);
            end
            
            xname = ['x', suffix, '_', num2str(obj.id)];
            mpol(xname, length(obj.vars.x), 1);
                        
            x_new = eval(xname);
            
            vars_new= struct('t', t_new, 'x', x_new);
            
            if nargin == 3
                vars_new.supp = subs_vars(supp_old, obj.get_vars(), ...
                                [t_new; x_new]);
            end
        end
        
        %% Constraints
        function cons = liou_con(obj, d, f)
            %LIOU_CON generate liouville constraint within location
            %
            %do not yet set this equal to zero (arithmetic operations not
            %are defined for constraints)
            if nargin == 2
                f = obj.f;
            end
            
            Ay_init = 0;
            if ~isempty(obj.meas_init)
                Ay_init =  obj.meas_init.mom_monom(d);
            end
            
            Ay_term = 0;
            if ~isempty(obj.meas_term)
                Ay_term = -obj.meas_term.mom_monom(d);
            end
            Ay_occ  =  obj.meas_occ.mom_lie(d, obj.vars, f);
            
            %TODO: Digital dynamics
            %Ay_occ  =  obj.meas_occ.mom_push(d, obj.vars, f);
            cons = Ay_init + Ay_term + Ay_occ;
        end
        
        function supp_con_out = supp_con(obj)
            
            if ~isempty(obj.meas_term)
                term_supp =  obj.meas_term.supp;
            else
                term_supp =  [];
            end
            
            if ~isempty(obj.meas_init)
                init_supp =  obj.meas_init.supp;
            else
                init_supp =  [];
            end
            
            supp_con_out = [init_supp;
                            term_supp;
                            obj.meas_occ.supp];
        end
        
        function [obj_max, obj_con] = objective_con(obj, d, objective)
            %OBJECTIVE_CON deal with the objective, which may be maximin
            if nargin == 2
                objective = obj.objective;
            end
                                    
            obj_con = [];
            if isempty(objective)
                obj_max = 0;
            elseif length(objective) == 1    
                obj_subs = obj.meas_term.var_sub(obj.vars, objective);
                obj_max = mom(obj_subs);                            
            else
                obj_subs = obj.meas_term.var_sub(obj.vars, objective);
                q_name = ['q_', num2str(obj.id)];
                mpol(q_name, 1, 1);
                q = eval(q_name);
                muq = meas(q);
                obj.cost_q = q;
                
                obj_max = mom(q);
                obj_con = [mass(q) == 1; (mom(q) <= mom(obj_subs));];
            end
            
        end
        
        function mass = mass_init(obj)
            mass = obj.meas_init.mass();
        end
        
        %% Recovery
        function s_out = mmat_corner(obj)
            s_out  = struct('init', [], 'term', [], 'occ', []);
            if ~isempty(obj.meas_init)
                s_out.init = obj.meas_init.mmat_corner();
            end
            if ~isempty(obj.meas_term)
                s_out.term = obj.meas_term.mmat_corner();
            end
            s_out.occ  = obj.meas_occ.mmat_corner();
        end
        
        function [optimal, mom_out, corner] = recover(obj, tol)
            %RECOVER if top corner of the moment matrix is rank-1, then
            %return approximate optimizer
            
            if nargin < 2
                tol = 5e-4;
            end
                        
            if isempty(obj.meas_init)
                opt_init = 1;
                mom_init.t = []; mom_init.x = [];
                corner_init = 0;
            else
                [opt_init, mom_init, corner_init] = obj.meas_init.recover(tol);
            end
            if isempty(obj.meas_term)
                opt_term = 1;
                mom_term.t = []; mom_term.x = [];
                corner_term = 0;
            else
                [opt_term, mom_term, corner_term] = obj.meas_term.recover(tol);
            end
            
            optimal = opt_init && opt_term;
            
            mom_out = struct('t0', mom_init.t, 'x0', mom_init.x, ...
                             'tp', mom_term.t, 'xp', mom_term.x);     
            corner = struct('init', corner_init, 'term', corner_term);
        end
        
        %% Dual variables
        
        function obj = dual_process(obj, order, v, beta, gamma)
             %DUAL_PROCESS turn the dual variables from solution into 
             %polynomials and interpretable quantities
             
             %numeric quantities
             obj.dual.solved = 1;
             
             obj.dual.beta = beta;
             obj.dual.gamma = gamma;
             
             %process the polynomials
             
             monom = mmon(obj.get_vars(), 0, 2*order);
             obj.dual.v = v'*monom;
             
             obj.dual.Lv = diff(obj.dual.v, obj.vars.x)*obj.f;
             if ~isempty(obj.vars.t)
                obj.dual.Lv = diff(obj.dual.v, obj.vars.t) + obj.dual.Lv;
             end
            
        end        
        
        %it may be worthwhile to set location_time_indep as its own class
        function f_out = f_eval(obj, t, x)
            %evaluate v
            if obj.TIME_INDEP
                f_out = eval(obj.f, obj.get_vars(), x);
            else 
                f_out = eval(obj.f, obj.get_vars(), [t; x]);
            end
        end
        
        function v_out = v_eval(obj, t, x)
            %evaluate v
            if obj.TIME_INDEP
                v_out = eval(obj.dual.v, obj.get_vars(), x);
            else
                v_out = eval(obj.dual.v, obj.get_vars(), [t; x]);
            end
        end
        
        function Lv_out = Lv_eval(obj, t, x)
            %evaluate Lv
            if obj.TIME_INDEP
                Lv_out = eval(obj.dual.Lv, obj.get_vars(), x);
            else
                Lv_out = eval(obj.dual.Lv, obj.get_vars(), [t; x]);
            end
        end
        
        function obj_out = obj_eval(obj, t, x)
            %evaluate objective
            if obj.TIME_INDEP
                obj_out = eval(obj.objective, obj.get_vars(), x);
            else
                obj_out = eval(obj.objective, obj.get_vars(), [t; x]);
            end
        end
        
        function supp_out = supp_eval(obj, t, x)
            %is (t, x) in the support of the location?
            if obj.TIME_INDEP
                supp_out =  all(eval(obj.supp, obj.get_vars(), [t; x]));
            else
                supp_out =  all(eval(obj.supp, obj.get_vars(), x));
            end
        end
        
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
        
        function cb = cost_beta(obj, t, x)
            if isempty(obj.objective)
                cb = zeros(size(t));
            else
                cb = eval(obj.dual.beta'*obj.objective, obj.get_vars(), [t; x]);
            end
        end

        function nn_out = nonneg(obj, t, x)
            %nonnegative functions at this location
            
            if isempty(obj.meas_init)
                nn_init = 0;
            else
                nn_init = obj.dual.gamma - obj.dual.v;
            end
            
            if isempty(obj.meas_term)
                nn_term = 0;
            else
                nn_term = obj.dual.v - obj.dual.beta'*obj.objective;
            end
            
            nn = [nn_init; -obj.dual.Lv; nn_term];
            
            if obj.TIME_INDEP
                nn_out = eval(nn, obj.get_vars(), x);                                   
            else
                nn_out = eval(nn, obj.get_vars(),  [t; x]);                                   
            end
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
        
        %% Sampler
        function out_sim = sample_traj_loc(obj, t0, x0, Tmax, curr_event)
            %SAMPLE_TRAJ_LOC Sample a single trajectory starting at (t0, x0) in
            %this location. Stop when the the trajectory hits a guard or
            %strays outside the location's support region
            %
            %curr_event handles the event detection for leaving the support
            %region, and guards if enabled.
            %
            %OUTPUT:
            %out_sim is a struct holding the simulation output: time,
            %state, objective, and nonnegative functions from the dual
            %solution of SDP.
            
            if nargin < 5
                curr_event = @obj.supp_event;
            end
            
            
            %simulate the trajectory
            curr_ode_options = odeset('Events',curr_event, 'RelTol', 1e-7, ...
                                      'AbsTol', 1e-8, 'MaxStep', 0.01);
        
            out_sim = struct;
            [out_sim.t, out_sim.x] = ode15s(@obj.f_eval, [t0, Tmax], x0, curr_ode_options);
            
            %evaluate nonnegative functions
            if obj.dual.solved
                out_sim.nonneg = obj.nonneg(out_sim.t', out_sim.x')';
            end
            
            out_sim.objective = obj.obj_eval(out_sim.t', out_sim.x')';
            out_sim.id = obj.id;
        end
                
    end
end

