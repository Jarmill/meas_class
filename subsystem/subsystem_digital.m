classdef subsystem_digital
    %SUBSYSTEM_DIGITAL A subsystem x+=f(t, x, th, w) of a possibly 
    %uncertain discrete dynamical system in measure analysis    
    
    properties
        vars = struct('t', [], 'x', [], 'th', [], 'w', []);
        
        %measures
        meas_occ = [];  %occupation measure          
        
        f = [];         %dynamics
           
        %variables and identification
        loc_id = [];
        sys_id = [];
        prefix = '';
        
        nn = [];        %vector of nonnegative functions
        
        dual = struct('v', 0, 'Lv', 0);
        
        supp;
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
            obj.supp = loc_supp;            
            obj.vars = loc_supp.vars;
            obj.supp.X_sys = obj.supp.get_X_sys_ind(sys_id);
            
            %occupation measure definition
            obj.f = f;
            
            obj.meas_occ  = obj.meas_def('occ');                                       
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
            %contribution of subsystem towards liouville constraints
            if obj.supp.DIGITAL
                Ay = obj.cons_push(d);
            else
                Ay = obj.cons_lie(d);
            end
        end
        
        %TODO: lie derivative (continuous)
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
        
        function Ay = cons_push(obj, d)
            %pushforward comparision affine combination
            %discrete systems only
            
            Ay = obj.meas_occ.mom_push(d, obj.vars, obj.f) - ...
                              obj.meas_occ.mom_monom_proj(d);
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
            obj.dual.Lv = diff(v, obj.vars.x)*obj.f;
            obj.dual.zeta = zeta;
            if ~isempty(obj.vars.t)
                obj.dual.Lv = obj.dual.Lv + diff(v, obj.vars.t);
            end
            
            %process the box dual variables           
            if isempty(zeta)
                obj.nn = -obj.dual.Lv;
            else
                Nb = length(obj.vars.b);
                %store all derivatives of v with respect to box occupation                
                Lv_box = zeros(Nb+1, 1)*obj.vars.b(1);
                
                for i = 1:(Nb+1)
                    Lv_box(i) = diff(v, obj.vars.x)*obj.f_box(i);
                    
                    if i==1 && ~isempty(obj.vars.t)
                        Lv_box(i) = Lv_box(i) + diff(v, obj.vars.t);
                    end
                end
                obj.Lv_box = Lv_box;
                %TODO: check the signs of these 
                nn_occ = -Lv_box(1) + sum(zeta);
                nn_box = -Lv_box(2:end) + zeta;
                nn_comp = zeta;
                
                obj.nn = [nn_occ; nn_box; nn_comp];
                
            end           
        end
        
        function obj = dual_process_digital(obj, v)
            %Dual process when dynamics are digital (no box)
            pushforward = eval(v, obj.vars.x, obj.f);
            Lv = pushforward - v;
            obj.nn = -Lv;
        end
        
        %% Evaluation and Sampling
        %TODO: Implement this
    end
end

