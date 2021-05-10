classdef subsystem_poly < subsystem_interface
    %SUBSYSTEM_POLY A subsystem x'=f(t, x, th, w, b) of a dynamical system 
    %where the uncertainty b is constrained to a polytope. This class is meant
    %may be applied to data-driven measure analysis
    
    properties
                
        %additional measures for polytopic uncertainty
        meas_pos  = {};  %box occupation measures
        meas_neg  = {};  %box occupation measures
        meas_comp = {};  %box-complement occupation measures                
        
        varnames = {'t', 'x', 'th', 'w'}; %names of variables in measure
        
        f_poly = [];     %affine decomposition of dynamics 
                        %{no input, input 1, input 2, ...}        
                        
        poly = [];
    end

    
    methods
        %% Constructor
        function obj = subsystem_poly(loc_supp, f, sys_id, loc_id)
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

                
                obj.poly = loc_supp.poly;
                
                Nconstraint = length(obj.poly.b);
                obj.meas_pos = cell(Nb, 1);
                obj.meas_neg = cell(Nb, 1);
                obj.meas_comp = cell(Nconstraint, 1);
                for i = 1:Nb
                    %positive and negative input measures
                    if Nb == 1
                        suffix_add = [];
                    else
                        suffix_add = ['_', num2str(i)];
                    end
                    obj.meas_pos{i}  = obj.meas_def({'t', 'x', 'th', 'w'}, ['_pos', suffix_add], obj.supp);
                    obj.meas_neg{i}  = obj.meas_def({'t', 'x', 'th', 'w'}, ['_neg', suffix_add], obj.supp);                                        
                end
                
                for i = 1:Nconstraint
                    %polytope complement measure
                    suffix_add = ['_', num2str(i)];
                    obj.meas_comp{i} = obj.meas_def({'t', 'x', 'th', 'w'}, ['_comp', suffix_add], obj.supp);
                    
                end    
                %process the dynamics f in terms of box dynamics
                %f_poly: [no input, input 1, input 2, ...]
                obj.f_poly = zeros(length(obj.f), Nb+1) * obj.vars.b(1);
                f0 = subs(obj.f, obj.vars.b, zeros(Nb, 1));                    
                obj.f_poly(:, 1) = f0;

                %each input channel at a time
                I = eye(Nb);

                for k = 1:Nb
                    obj.f_poly(:, k+1) = subs(obj.f, obj.vars.b, I(:, k)) - f0;                        
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
                %base occupation measure (with no polytopic disturbances)
                Ay = obj.meas_occ.mom_lie(d, obj.get_vars, obj.f_poly(:, 1));
                
                
                for k = 1:Nb
                    Ay_pos = obj.meas_pos{k}.mom_lie(d, obj.get_vars, obj.f_poly(:, k+1), 0);
                    Ay_neg= obj.meas_neg{k}.mom_lie(d, obj.get_vars, obj.f_poly(:, k+1), 0);

                    %add contribution to lie derivative
                    Ay = Ay + Ay_pos - Ay_neg;
                end
                
            end
            
        end       
        
        function Ay = abscont_box(obj, d)
            %ABSCONT_BOX absolute continuity constraints of each input+complement with 
            %respect to the occupation measure
            Ay = [];
            
            Nconstraints = length(obj.poly.b);                        
            
            %cons b >= ans
            %-cons b <= -ans
            
            cons_mat = obj.poly.A;
            ans_vec = obj.poly.b;
            
            %moments of each measure
            mom_occ = obj.meas_occ.mom_monom(d);
            Nb = length(obj.vars.b);
            
            for j = 1:Nconstraints        
                
                %For each constraint Ab <= d
                %the abscont constraint is
                %Sum_i <v A_ij, (bi+)+(bi-)>  + <v, (slackj) = d_j <v, occ>
                
                poly_contrib = 0;
                for i = 1:Nb
                    weight_curr = cons_mat(j, i);
                    if weight_curr ~= 0

                        mom_pos= obj.meas_pos{i}.mom_monom(d) * weight_curr;
                        mom_neg= obj.meas_neg{i}.mom_monom(d) * weight_curr;
                    
                        poly_contrib = poly_contrib + mom_pos + mom_neg;
                    end
                end

                mom_comp = obj.meas_comp{j}.mom_monom(d);
                
                %absolute continuity constraint
                %A(input) <= b
                Ay_curr = -mom_occ * ans_vec(j) + poly_contrib + mom_comp;
                Ay = [Ay; Ay_curr]; 
            end
        end        
        

              
        %% Getters
        function supp_all = get_supp(obj)
            %SUPP_ALL: get support set of all measures in subsystem
            
            supp_all = obj.meas_occ.supp;
            
            %input measures
            for i = 1:length(obj.meas_pos)
                supp_pos = obj.meas_pos{i}.supp;
                supp_neg = obj.meas_neg{i}.supp;
                                
                supp_all = [supp_all; supp_pos; supp_neg];                
            end
            
            %complement measures
            for i = 1:length(obj.meas_comp)
                supp_comp = obj.meas_comp{i}.supp;
                supp_all = [supp_all; supp_comp];
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
                obj.dual.nn = -obj.dual.Lv;
            else
                Nb = length(obj.vars.b);
                Nconstraints = length(obj.poly.b);
                %store all derivatives of v with respect to box occupation                
                Lv_box = zeros(Nb+1, 1)*obj.vars.b(1);
                
                %TODO: fix this
                
                for i = 1:(Nb+1)
                    Lv_box(i) = diff(v, obj.vars.x)*obj.f_poly(:, i);
                    
                    if i==1 && ~isempty(obj.vars.t)
                        Lv_box(i) = Lv_box(i) + diff(v, obj.vars.t);
                    end
                end
                obj.dual.Lv_box = Lv_box;
                %TODO: check the signs of these 
                nn_occ = -Lv_box(1) + obj.poly.b'*(zeta);
                nn_poly = -Lv_box(2:end) + obj.poly.A'*zeta;
                nn_comp = zeta;
                
                obj.dual.nn = [nn_occ; nn_poly; nn_comp];                
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

