classdef subsystem < subsystem_interface
    %SUBSYSTEM A subsystem x'=f(t, x, th, w, b) of a possibly uncertain 
    %dynamical system in measure analysis    
    
    properties
                
        %additional measures for box uncertainty
        meas_box = {};  %box occupation measures
        meas_comp = {}; %box-complement occupation measures                
        
        varnames = {'t', 'x', 'th', 'w'}; %names of variables in measure
        
        f_box = [];     %affine decomposition of dynamics 
                        %{no input, input 1, input 2, ...}     
                        
    end

    
    methods
        %% Constructor
        function obj = subsystem(loc_supp, f, sys_id, loc_id)
            %SUBSYSTEM Construct a continuous (possibly uncertain) 
            %subsystem, fill in information
            
            %process input
            if nargin < 3
                sys_id = 1;
            end
            
            if nargin < 4
                loc_id = [];
            end
            
            %superclass constructor
            obj@subsystem_interface(loc_supp, f, sys_id, loc_id, @meas_uncertain);
%             obj.meas_type = @meas_uncertain;
            
            obj.dual = struct('v', 0, 'Lv', 0, 'Lv_box', 0, 'zeta', 0, 'nn', 0);
                       
            %box-occupation measures definition
            if ~isempty(obj.vars.b)
                Nb = length(obj.vars.b);
                obj.meas_box = cell(Nb, 1);
                obj.meas_comp = cell(Nb, 1);
                for i = 1:Nb
                    %box measure
                    if Nb == 1
                        suffix_add = [];
                    else
                        suffix_add = ['_', num2str(i)];
                    end
                    obj.meas_box{i}  = obj.meas_def({'t', 'x', 'th', 'w'}, ['_box', suffix_add], obj.supp);
                    
                    
                    %box complement measure
                    obj.meas_comp{i} = obj.meas_def({'t', 'x', 'th', 'w'}, ['_comp', suffix_add], obj.supp);
                    
                end    
                %process the dynamics f in terms of box dynamics
                %f_box: [no input, input 1, input 2, ...]
                obj.f_box = zeros(length(obj.f), Nb+1) * obj.vars.b(1);
                f0 = subs(obj.f, obj.vars.b, zeros(Nb, 1));                    
                obj.f_box(:, 1) = f0;

                %each input channel at a time
                I = eye(Nb);

                for k = 1:Nb
                    obj.f_box(:, k+1) = subs(obj.f, obj.vars.b, I(:, k)) - f0;                        
                end
                                                                
                
            end            
        end        
        
        %% Getters
        
        function vars_out = get_vars(obj)
            %GET_VARS_BOX include box variables b
            vars_out = [obj.vars.t; obj.vars.x; obj.vars.th; obj.vars.w];
        end
        
        function vars_out = get_vars_box(obj)
            %GET_VARS_BOX include box variables b
            vars_out = [obj.vars.t; obj.vars.x; obj.vars.th; obj.vars.w; obj.vars.b];
        end
        
        %% Constraints        
        
       
        function Ay = cons_liou(obj, d)
            %CONS_LIOU Liouville Equation includes an affine combination of
            %Lie derivatives (continuous systems only)
            
            if isempty(obj.vars.b)
                %no box inputs, simple to perform
                 Ay = obj.meas_occ.mom_lie(d, obj.get_vars, obj.f);
            else
                %non-trivial box inputs, more involved processing
                
                Nb = length(obj.vars.b);
                %base occupation measure (with no box disturbance)
%                 vars_
                Ay = obj.meas_occ.mom_lie(d, obj.get_vars, obj.f_box(:, 1));
                
                %each input channel at a time
%                 I = eye(Nb);
                
                for k = 1:Nb
%                     fk = subs(obj.f, obj.vars.b, I(:, k)) - f0;
                    Ay_curr = obj.meas_box{k}.mom_lie(d, obj.get_vars, obj.f_box(:, k+1), 0);
                    
                    %add contribution to lie derivative
                    Ay = Ay + Ay_curr;
                end
                
            end
            
        end       
        
        function Ay = abscont_box(obj, d)
            %ABSCONT_BOX absolute continuity constraints of each box+complement with 
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
                obj.dual.Lv = obj.Tmax*obj.dual.Lv + diff(v, obj.vars.t);
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
                        Lv_box(i) = obj.Tmax*Lv_box(i) + diff(v, obj.vars.t);
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

    end
end

