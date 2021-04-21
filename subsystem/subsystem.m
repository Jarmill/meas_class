classdef subsystem
    %SUBSYSTEM A subsystem x'=f(t, x, th, w, b) of a possibly uncertain 
    %dynamical system in measure analysis    
    
    properties
        vars = struct('t', [], 'x', [], 'th', [], 'w', [], 'b', []);
        
        %measures
        meas_occ = [];  %occupation measure
        meas_box = {};  %box occupation measures
        meas_comp = {}; %box-complement occupation measures                
        
        f = [];         %dynamics
        f_box = {};     %affine decomposition of dynamics 
                        %{no input, input 1, input 2, ...}
        
        %variables and identification
        loc_id = [];
        sys_id = [];
        prefix = '';                
        
        dual = struct('v', 0, 'Lv', 0, 'Lv_box', 0, 'zeta', 0, 'nn', 0);
        
        supp;       %support of measures
        supp_sys;   %support of only x
        DIGITAL = 0;
    end
    
    properties(Access = private)
        %use private properties for numeric function evaluations
        %object-oriented matlab is slow with public variables 
        f_;
        supp_sys_ = []; %support of state
        nn_ = 0;        %nonnegative function
    end
    
    methods
        %% Constructor
        function obj = subsystem(loc_supp, f, sys_id, loc_id)
            %SUBSYSTEM Construct a subsystem, fill in information
            
            %process input
            if nargin < 3
                sys_id = 1;
            end
            
            if nargin < 4
                loc_id = [];
            end
            
            %identification and names
            obj.loc_id = loc_id;
            obj.sys_id = sys_id;
            obj.prefix = ['_', num2str(sys_id), '_'];
            if ~isempty(loc_id)
                obj.prefix = [num2str(loc_id), obj.prefix];
            end
            
            %support            

            obj.vars = loc_supp.vars;
            obj.supp_sys= loc_supp.get_X_sys_ind(sys_id);
            obj.supp_sys_ = obj.supp_sys;
            obj.supp = loc_supp.supp_sys_pack(obj.supp_sys);           
            %occupation measure definition
            obj.f = f;
            obj.f_ = f; 
            
            obj.meas_occ  = obj.meas_def('occ');          
            
            %box-occupation measures definition
            if ~isempty(obj.vars.b)
                Nb = length(obj.vars.b);
                obj.meas_box = cell(Nb, 1);
                obj.meas_comp = cell(Nb, 1);
                for i = 1:Nb
                    %box measure
                    obj.meas_box{i}  = obj.meas_def(['box_', num2str(i)]);
                    
                    
                    %box complement measure
                    obj.meas_comp{i} = obj.meas_def(['comp_', num2str(i)]);
                    
                    
                    %process the dynamics f in terms of box dynamics
                    %f_box: {no input, input 1, input 2, ...}
%                     obj.f_box = cell(Nb+1, 1);
                    obj.f_box = zeros(Nb+1, 1) * obj.vars.b(1);
                    f0 = subs(obj.f, obj.vars.b, zeros(Nb, 1));                    
                    obj.f_box(:, 1) = f0;
                    
                    %each input channel at a time
                    I = eye(Nb);

                    for k = 1:Nb
                        obj.f_box(:, k+1) = subs(obj.f, obj.vars.b, I(:, k)) - obj.f_box(1);                        
                    end
                end                                                
                
            end            
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
                
                supp_new = subs_vars(obj.supp, obj.get_vars(), ...
                                [vars_new.t; vars_new.x; vars_new.th; vars_new.w]);
           
            
            %define the measure
            meas_new = meas_uncertain(vars_new, supp_new);
        end
        
        %% Getters
        
        function vars_out = get_vars(obj)
            %GET_VARS add more variables as necessary
            vars_out = [obj.vars.t; obj.vars.x; obj.vars.th; obj.vars.w];
        end
        
        function vars_out = get_vars_box(obj)
            %GET_VARS_BOX include box variables b
            vars_out = [obj.vars.t; obj.vars.x; obj.vars.th; obj.vars.w; obj.vars.b];
        end
        
        %% Constraints        
        function Ay = cons_liou(obj, d)
            %contribution of subsystem towards liouville constraints
            if obj.DIGITAL
                Ay = obj.cons_push(d);
            else
                Ay = obj.cons_lie(d);
            end
        end
       
        function Ay = cons_lie(obj, d)
            %lie derivative affine combination
            %continuous systems only
            
%             Ay = 0;
            if isempty(obj.vars.b)
                %no box inputs, simple to perform
                 Ay = obj.meas_occ.mom_lie(d, obj.vars, obj.f);
            else
                %non-trivial box inputs, more involved processing
                
                Nb = length(obj.vars.b);
                %base occupation measure (with no box disturbance)
                %TODO: replace with f_box
%                 f0 = subs(obj.f, obj.vars.b, zeros(Nb, 1));
                Ay = obj.meas_occ.mom_lie(d, obj.vars, obj.f_box(:, 1));
                
                %each input channel at a time
%                 I = eye(Nb);
                
                for k = 1:Nb
%                     fk = subs(obj.f, obj.vars.b, I(:, k)) - f0;
                    Ay_curr = obj.meas_box{k}.mom_lie(d, obj.vars, obj.f_box(:, k+1), 0);
                    
                    %add contribution to lie derivative
                    Ay = Ay + Ay_curr;
                end
                
            end
            
        end       
        
        function Ay = abscont(obj, d)
            %absolute continuity constraints of each box+complement with 
            %respect to the occupation measure
            Ay = [];
            
            %moments of each measure
            mom_occ = obj.meas_occ.mom_monom(d);
            for i = 1:length(obj.vars.b)
                mom_box  = obj.meas_box{i}.mom_monom(d);
                mom_comp = obj.meas_comp{i}.mom_monom(d);
                
                %absolute continuity constraint
                Ay_curr = -mom_occ + mom_box + mom_comp;
                Ay = [Ay; Ay_curr]; 
            end
        end
        
        function Ay = cons_push(obj, d)
            %pushforward comparision affine combination
            %discrete systems only
            
            Ay = obj.meas_occ.mom_push(d, obj.vars, obj.f) - ...
                              obj.meas_occ.mom_monom_proj(d);
        end
        

              
        %% Getters
        function supp_all = get_supp(obj)
            %SUPP_ALL: get support set of all measures in subsystem
            
            supp_all = obj.meas_occ.supp;
            for i = 1:length(obj.meas_box)
                supp_box = obj.meas_box{i}.supp;
                supp_comp = obj.meas_comp{i}.supp;
                
                supp_all = [supp_all; supp_box; supp_comp];
            end
        end
        
        
        %% Dual Recovery 
        function obj = dual_process(obj, v, zeta)
            %DUAL_PROCESS store dual functions and compute nonnegative
            %functions for this subsystem
            
            %auxiliary function v            
            obj.dual.v = v;
            
            if obj.DIGITAL
                %Dual process when dynamics are digital (no box)
                pushforward = eval(v, obj.vars.x, obj.f);
                Lv = pushforward - v;
            else
                obj.dual.Lv = diff(v, obj.vars.x)*obj.f;
                obj.dual.zeta = zeta;
                if ~isempty(obj.vars.t)
                    obj.dual.Lv = obj.dual.Lv + diff(v, obj.vars.t);
                end
            end
            
            %process the box dual variables           
            if isempty(zeta)
                obj.dual.nn = -obj.dual.Lv;
            else
                Nb = length(obj.vars.b);
                %store all derivatives of v with respect to box occupation                
                Lv_box = zeros(Nb+1, 1)*obj.vars.b(1);
                
                for i = 1:(Nb+1)
                    Lv_box(i) = diff(v, obj.vars.x)*obj.f_box(:, i);
                    
                    if i==1 && ~isempty(obj.vars.t)
                        Lv_box(i) = Lv_box(i) + diff(v, obj.vars.t);
                    end
                end
                obj.dual.Lv_box = Lv_box;
                %TODO: check the signs of these 
                nn_occ = -Lv_box(1) + sum(zeta);
                nn_box = -Lv_box(2:end) + zeta;
                nn_comp = zeta;
                
                obj.dual.nn = [nn_occ; nn_box; nn_comp];                
            end          
            
            obj.nn_ = obj.dual.nn;
        end
        

        
        %% Function Evaluation (for sampling)
        %TODO: Implement this
        %nonnegativity evaluation: include supports        
        
        %A function to evaluate dynamics f_
        function f_out = f_eval(obj, data)
            %data: [t, x, th, w, b] as required            
            f_out = eval(obj.f_, obj.get_vars_box(), data);
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
end

