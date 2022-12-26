classdef location_interface < handle
    %LOCATION_INTERFACE A location (space) of a dynamical system
    %   includes descriptions of the space as well as measures
    
    properties
        %location number
        id = [];
        
        %measures
        init;   %initial measure
        term;   %terminal measure
        sys;    %subsystems (occupation measures)
        
        %variables
        vars;
        cost_q = []; %multiple costs
        supp;   %support
        
        f;          %dynamics
        objective;
        
        %dual variable recovery
        len_dual;   %length of coefficients of each polynomial
        dual;       %dual polynomials, including nonnegative functions
        %dual = struct('v', 0, 'beta', 1, 'gamma', 0, 'nn', 0, 'solved', 0);         
                
        %bound by loc_sampler constructor
        sampler = [];
        
    end
    
    properties(Access = protected)
        %use private properties for numeric function evaluations
        %object-oriented matlab is slow with public variables 
        %
        %TODO: check if protected is also slow
        %this is copied from subsystem/subsystem_interface
        TIME_INDEP = 0;
        supp_X_;
        nn_;
    end
    
    methods
        function obj = location_interface(loc_supp, f, objective, id)
            %LOCATION_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.supp = loc_supp;
            obj.vars = loc_supp.vars;
            
            obj.objective = objective;
            obj.id = id;
            
            %dynamics
            if ~iscell(f)
                obj.f = {f};
            else
                obj.f = f;            
            end
            
            %systems (occupation measures)
            Nsys = length(obj.f);

            
            if ~obj.supp.TIME_INDEP % && obj.supp.SCALE_TIME
                %scale for time if time is a variable
                Tmax = obj.supp.Tmax;
                for i = 1:Nsys
                    obj.f{i} = obj.f{i}*Tmax;
                end
                obj.supp.Tmax = 1;
            end
            
            %supports
            if ~iscell(obj.supp.X_sys)
                obj.supp.X_sys = {obj.supp.X_sys};         
            end
            
            obj.supp_X_ = obj.supp.X;
            
            %initial measures
            if ~isempty(obj.supp.X_init) || ~isempty(obj.supp.mom_init)
%                 obj.init = obj.meas_def_end('0', obj.supp.supp_init());
                obj.init = meas_init(obj.supp, obj.id);
            end            
            
            %terminal measures
            if ~isempty(objective)
%                  obj.term = obj.meas_def_end('p', obj.supp.supp_term());                
                obj.term = meas_term(obj.supp, obj.id);             
            end       
            
            %occupation measures get created in the subclass location or
            %similar subclass
        end     
        
        
        %% getters
        function vars_out = var_index(obj, vars_in, varnames)
                   
            %generate a stack of variables
            %inefficient code, fix later
            vars_out = struct;
            for i = 1:length(varnames)
                curr_var = varnames{i};
                if isempty(vars_in.(curr_var))
                    vars_out.(curr_var) = [];
                else
                    vars_out.(curr_var) = reshape(vars_in.(curr_var), [], 1);
                end
            end
        end
        
        %% moments
        function mass_out = mass_occ(obj)
            %return the sum of the masses of all occupation measures
            %useful for constraining a time-independent free-time system to
            %have finite time
            mass_out = 0;
            for i =1:length(obj.sys)
                mass_out = mass_out + obj.sys{i}.mass_occ();
            end
        end
        
        function mass = mass_init(obj)
            mass = obj.init.mass();
        end
        
        %% support
        function supp_con_out = supp_con(obj)
            %SUPP_CON get support constraints of measures
            
            
            %terminal measure support 
            if ~isempty(obj.term)
                term_supp =  obj.term.supp();
            else
                term_supp =  [];
            end
            
            %initial measure support 
            if ~isempty(obj.init)
                init_supp =  obj.init.supp();
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
        
        %% constraints
        function cons = liou_con(obj, d)
            %LIOU_CON generate liouville constraint within location
            %
            %do not yet set this equal to zero (arithmetic operations not
            %are defined for constraints)
            
            Ay_init = 0;
            if ~isempty(obj.init)
                Ay_init =  obj.init.mom_monom(d);
            end
            
            Ay_term = 0;
            if ~isempty(obj.term)
                Ay_term = -obj.term.mom_monom(d);
            end
            
            %TODO replace with a subsystem call
            Ay_occ = 0;
            for i = 1:length(obj.sys)

                Ay_curr = obj.sys{i}.cons_liou(d);
                Ay_occ  =  Ay_occ + Ay_curr;
            end
            
            cons = Ay_init + Ay_term + Ay_occ;
        end
        
        function [cons, len_abscont] = abscont_box_con(obj, d)
            %ABSCONT_BOX_CON constraint for absolute continuity in box-disturbance
            cons = [];
            len_abscont = zeros(length(obj.sys), 1);
            for i = 1:length(obj.sys)
                curr_abscont = obj.sys{i}.abscont_box(d);
                cons = [cons; curr_abscont];
                len_abscont(i) = length(curr_abscont);
            end
        end
        
        %% Recovery
        function s_out = mmat_corner(obj)
            s_out  = struct('init', [], 'term', [], 'occ', []);
            if ~isempty(obj.init)
                s_out.init = obj.init.mmat_corner();
            end
            if ~isempty(obj.term)
                s_out.term = obj.term.mmat_corner();
            end
%             s_out.occ  = obj.meas_occ.mmat_corner();
        end
        
        function [optimal, mom_out, corner] = recover(obj, tol)
            %RECOVER if top corner of the moment matrix is rank-1, then
            %return approximate optimizer
            
            if nargin < 2
                tol = 5e-4;
            end
                        
            if isempty(obj.init)
                opt_init = 1;
                mom_init.t = []; mom_init.x = [];
                corner_init = 0;
            else
                [opt_init, mom_init, corner_init] = obj.init.recover(tol);
            end
            if isempty(obj.term)
                opt_term = 1;
                mom_term.t = []; mom_term.x = [];
                corner_term = 0;
            else
                [opt_term, mom_term, corner_term] = obj.term.recover(tol);
            end
            
            optimal = opt_init && opt_term;
            
            mom_out = struct('t0', mom_init.t, 'x0', mom_init.x, ...
                             'tp', mom_term.t, 'xp', mom_term.x);     
            corner = struct('init', corner_init, 'term', corner_term);
        end
        
        %% Sampling
        function obj_out = obj_eval(obj, t, x)
            %evaluate objective
            if obj.TIME_INDEP
                obj_out = eval(obj.objective, obj.vars.x, x);
            else
                obj_out = eval(obj.objective, [obj.vars.t; obj.vars.x], [t'; x']);
            end
        end
        
        function supp_out = supp_eval(obj, t, x)
            %is (t, x) in the support of the location?
%             if obj.loc_supp.TIME_INDEP
%                 supp_out =  all(eval(obj.supp.X, obj.get_vars(), [t; x]));
%             else
                supp_out =  all(eval(obj.supp_X_, obj.vars.x, x));
%             end
        end

        function x_proj = supp_proj(obj, x)
            %project a point x onto the support set
            %used to avoid/throw under rug numerical difficulties in
            %sampling routines

            %should be overloaded in child classes
            x_proj = x;
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
        
        %% overloads
        function e = isempty(obj)
            %is the support empty?
            %as in supp = []. The harder question would be 'does the basic
            %semialgebraic set formed by the constraints satisfy a
            %nullstellensatz?'
            e = isempty(obj.supp);
        end
        
    end
    
    methods(Abstract)
        all_cons(obj, d)
        %ALL_CONS all constraints involving solely this location        
        
        
        [len_out] = len_eq_cons(obj)
        %LEN_EQ_CONS Number of equality constraints strictly in this
        %location 
        
        objective_con(obj, objective)
        %OBJECTIVE_CON deal with the objective, which may be maximin
        
        dual_process(obj, d, rec_eq, rec_ineq, gamma)
         %DUAL_PROCESS turn the dual variables from solution into 
         %polynomials and interpretable quantities
    end
end

