classdef location < handle
    %LOCATION A location (space) of a hybrid system
    %   includes descriptions of the space as well as measures
    
    properties
        %location number
        id; 
        
        %measures
        meas_init;
        meas_term;
        sys;
        
        %variables
%         vars = struct('t', [], 'x', []);
        vars = struct('t', [], 'x', [], 'th', [], 'w', [], 'b', []);
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
        function obj = location(id, loc_supp, f, objective)
            %Location Construct an instance of this class
            %   Detailed explanation goes here
            
            %fill in properties
            obj.id = id;
%             obj.vars = vars;
%             if ~isfield(vars, 't') || isempty(vars.t)
%                 obj.vars.t = []; %time-independent
%                 obj.loc_supp.TIME_INDEP = 1;            
%             end

            
            obj.supp = loc_supp;
            obj.vars = loc_supp.vars;
            
            %dynamics
            if ~iscell(obj.f)
                obj.f = {f};
            else
                obj.f = f;            
            end
           
            
            if ~iscell(obj.supp.X_sys)
                obj.supp.X_sys = {obj.supp.X_sys};
            else
                obj.supp.X_sys = obj.supp.X_sys;            
            end
                       
            if nargin < 4                  
                %by default, no objective
                obj.objective = [];
            else
                obj.objective = objective;
            end
            

            %systems (occupation measures)
            Nsys = length(obj.f);

            
            if ~obj.supp.TIME_INDEP
                %scale for time if time is a variable
                Tmax = obj.supp.Tmax;
                for i = 1:Nsys
                    obj.f{i} = obj.f{i}*Tmax;
                end
                obj.supp.Tmax = 1;
            end
            
            obj.sys = cell(Nsys, 1);
            %subsystems
            for i = 1:Nsys
                obj.sys{i} = subsystem(obj.supp, obj.f{i}, i, id);
            end
            
            
            %TODO: add facility for multiple initial/final regions (unions)
            
            %initial measures
            if ~isempty(obj.supp.X_init)
                obj.meas_init = obj.meas_def_end('0', obj.supp.supp_init());
            end
            
            %terminal measures
            if ~isempty(objective)
                 obj.meas_term = obj.meas_def_end('p', obj.supp.supp_term());                
            end
            

                       
        end
        
        function vars_out = get_vars_end(obj)
            %GET_VARS_END variables at endpoint measures
            %   initial and terminal, without time-dependent
            vars_out = [obj.vars.t; obj.vars.x; obj.vars.th];
        end
        
        function vars_out = get_vars(obj)
            %GET_VARS add more variables as necessary
            vars_out = [obj.vars.t; obj.vars.x; obj.vars.th; obj.vars.w];
        end
        
        %% measure definition
        function meas_new = meas_def_end(obj, suffix, supp_ref)           
            %declare a variable for each initial or terminal measure
            %does not include 'w' or 'b' (time-dependent inputs)
            vars_new = struct('t', [], 'x', [], 'th', []);           
            varnames = fields(vars_new);
            for i = 1:length(varnames)
                curr_name = varnames{i};
                curr_var = obj.vars.(curr_name);
                
                if ~isempty(curr_var)
                    %declare a new variable
                    new_name = [curr_name, suffix];
                    mpol(new_name, length(curr_var), 1);
                    %load the new variable into vars_new
                    vars_new.(curr_name) = eval(new_name);
                end
%                 obj.vars.(curr_var) = vars.(curr_var);
            end
            
           
                %create new support as well
%                 supp_ref = obj.supp.supp_sys_pack(obj.supp.X_sys);
                supp_new = subs_vars(supp_ref, [obj.vars.t; obj.vars.x; obj.vars.th], ...
                                [vars_new.t; vars_new.x; vars_new.th]);
           
            
            %define the measure
            meas_new = meas_base(vars_new, supp_new);
        end
        
        
        %% Constraints
        function cons = liou_con(obj, d)
            %LIOU_CON generate liouville constraint within location
            %
            %do not yet set this equal to zero (arithmetic operations not
            %are defined for constraints)
            
            Ay_init = 0;
            if ~isempty(obj.meas_init)
                Ay_init =  obj.meas_init.mom_monom(d);
            end
            
            Ay_term = 0;
            if ~isempty(obj.meas_term)
                Ay_term = -obj.meas_term.mom_monom(d);
            end
            
            %TODO replace with a subsystem call
            Ay_occ = 0;
            for i = 1:length(obj.sys)

                Ay_curr = obj.sys{i}.cons_liou(d);
                Ay_occ  =  Ay_occ + Ay_curr;
            end
            
            cons = Ay_init + Ay_term + Ay_occ;
        end
        
        function [cons, len_abscont] = abscont_con(obj, d)
            %constraint for absolute continuity in box-disturbance
            cons = [];
            len_abscont = zeros(length(obj.sys), 1);
            for i = 1:length(obj.sys)
                curr_abscont = obj.sys{i}.abscont(d);
                cons = [cons; curr_abscont];
                len_abscont(i) = length(curr_abscont);
            end
        end
        
        
        function supp_con_out = supp_con(obj)
            %SUPP_CON get support constraints of measures
            
            
            %terminal measure support 
            if ~isempty(obj.meas_term)
                term_supp =  obj.meas_term.supp;
            else
                term_supp =  [];
            end
            
            %initial measure support 
            if ~isempty(obj.meas_init)
                init_supp =  obj.meas_init.supp;
            else
                init_supp =  [];
            end
            
            %subsystem measure support
            sys_supp = [];
            for i = 1:length(obj.sys)
                sys_supp = [sys_supp; obj.sys{i}.get_supp()];
            end
            
            supp_con_out = [init_supp;
                            term_supp;
                            sys_supp];
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
                obj_subs = obj.meas_term.var_sub_end(obj.vars, objective);
                obj_max = mom(obj_subs);                            
            else
                obj_subs = obj.meas_term.var_sub_end(obj.vars, objective);
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
%             s_out.occ  = obj.meas_occ.mmat_corner();
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
             
             %TODO: fix this so that boxes can be used
             
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
            if obj.loc_supp.TIME_INDEP
                f_out = eval(obj.f, obj.get_vars(), x);
            else 
                f_out = eval(obj.f, obj.get_vars(), [t; x]);
            end
        end
        
        function v_out = v_eval(obj, t, x)
            %evaluate v
            if obj.loc_supp.TIME_INDEP
                v_out = eval(obj.dual.v, obj.get_vars(), x);
            else
                v_out = eval(obj.dual.v, obj.get_vars(), [t; x]);
            end
        end
        
        function Lv_out = Lv_eval(obj, t, x)
            %evaluate Lv
            if obj.loc_supp.TIME_INDEP
                Lv_out = eval(obj.dual.Lv, obj.get_vars(), x);
            else
                Lv_out = eval(obj.dual.Lv, obj.get_vars(), [t; x]);
            end
        end
        
        function obj_out = obj_eval(obj, t, x)
            %evaluate objective
            if obj.loc_supp.TIME_INDEP
                obj_out = eval(obj.objective, obj.get_vars(), x);
            else
                obj_out = eval(obj.objective, obj.get_vars(), [t; x]);
            end
        end
        
        function supp_out = supp_eval(obj, t, x)
            %is (t, x) in the support of the location?
%             if obj.loc_supp.TIME_INDEP
%                 supp_out =  all(eval(obj.supp.X, obj.get_vars(), [t; x]));
%             else
                supp_out =  all(eval(obj.supp.X, obj.vars.x, x));
%             end
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
            
            if obj.loc_supp.TIME_INDEP
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

