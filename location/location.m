classdef location < location_interface
    %LOCATION A location (space) of a dynamical system
    %   includes descriptions of the space as well as measures
    %   used for continuous or discrete time ODE systems
    
    properties
        %variables
%         vars = struct('t', [], 'x', []);
%         vars = struct('t', [], 'x', [], 'th', [], 'w', [], 'b', []);
        varnames = {'t','x','th','w','b'};
        

        
           
    end
       
    methods
        function obj = location(loc_supp, f, objective, id)
            %Location Construct an instance of this class
            %   Detailed explanation goes here
            
            %fill in properties
            
            if nargin < 3             
                %by default, no objective
                objective = [];                
            end
            
            if nargin < 4
                id = [];            
            end
            obj@location_interface(loc_supp, f, objective, id);
            
            %TODO: make id the last argument                                       
                       
            


            Nsys = length(obj.f);
            obj.sys = cell(Nsys, 1);
            %subsystems
            for i = 1:Nsys                
                if obj.supp.DIGITAL
                    obj.sys{i} = subsystem_digital(obj.supp, obj.f{i}, i, id);
                else
                    obj.sys{i} = subsystem(obj.supp, obj.f{i}, i, id);
                end
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
        
        function vars_out = get_vars_box(obj)
            %GET_VARS_BOX include box variables b
            vars_out = [obj.vars.t; obj.vars.x; obj.vars.th; obj.vars.w; obj.vars.b];
        end
        
        
        %% Constraints
        
        function [objective, cons_eq, cons_ineq, len_dual] = all_cons(obj, d)
            %ALL_CONS all constraints involving solely location
            %does not include sum of mass of initial measures
            %Output:
            %   cons_eq: equality constraints
            %   cons_ineq: inequality constraints (objective)
            
            %gather all constraints together
            liou = obj.liou_con(d);
            len_liou = length(liou);
            [abscont_box, len_abscont] = obj.abscont_box_con(d);
            
            [objective, cons_ineq] = obj.objective_con();
            
            %package up the output
            len_dual = struct;
            len_dual.v = len_liou;
            len_dual.zeta = len_abscont;
            len_dual.beta = length(cons_ineq);
            
            %ensure this iss the correct sign
            cons_eq = [-liou; abscont_box]==0;                        
        end                              
       
        
        function [len_out] = len_eq_cons(obj)
            %LEN_EQ_CONS Number of equality constraints strictly in this
            %location 
            len_out = obj.len_dual.v + sum(obj.len_dual.zeta);
        end
        
        function [obj_max, obj_con] = objective_con(obj, objective)
            %OBJECTIVE_CON deal with the objective, which may be maximin
                                    
            %TODO: This should maybe go in the manager
            %The current implementation is only for peak estimation
            
            %TODO: include support for putting objectives on initial and
            %occupation measures as well as the terminal measure
            if nargin == 1
                objective = obj.objective;
            end
                                    
            obj_con = [];
            
            var_end = obj.var_index(obj.vars, {'t', 'x', 'th'});
            if isempty(objective)
                obj_max = 0;
            elseif length(objective) == 1    
                obj_subs = obj.term.var_sub_mom(var_end, objective);
                obj_max = (obj_subs);                            
            else
                obj_subs = obj.term.var_sub_mom(var_end, objective);
                q_name = ['q_', num2str(obj.id)];
                mpol(q_name, 1, 1);
                q = eval(q_name);
                muq = meas(q);
                obj.cost_q = q;
                
                obj_max = mom(q);
                obj_con = [mass(q) == 1; (mom(q) <= mom(obj_subs));];
            end            
        end
               
        
        %% Dual variables
        
        function obj = dual_process(obj, d, rec_eq, rec_ineq, gamma)
             %DUAL_PROCESS turn the dual variables from solution into 
             %polynomials and interpretable quantities
             %
             %Input:
             %  d:          2*order, order of auxiliary polynomials
             %  rec_eq:     dual variables from equality constraints
             %  rec_ineq:   dual variables from inequality constraints
             %  gamma:      objective value (as a dual variable)
             %TODO: fix this so that boxes can be used
             
             %numeric quantities
             obj.dual.solved = 1;
             
             obj.dual.beta = rec_ineq;
             obj.dual.gamma = gamma;
             
             %process the polynomials
             
             %auxiliary function v
             v_coeff = rec_eq(1:obj.len_dual.v);
             monom = mmon(obj.get_vars_end(), 0, d);
             obj.dual.v = v_coeff'*monom;
             
             count_zeta = obj.len_dual.v;
             
             %iterate through all subsystems
             Nb = length(obj.vars.b);
             monom_all = mmon(obj.get_vars(), 0, d);
             
             %TODO: confirm that all abscont relations have the same length
             len_monom_all = length(monom_all);
             for i = 1:length(obj.sys)       
                 
                 %untangle the box zeta functions
                 zeta = [];
                 for j = 1:Nb
                     zeta_coeff = rec_eq(count_zeta + (1:len_monom_all));
                     zeta_curr = zeta_coeff'*monom_all;
                     zeta = [zeta; zeta_curr];   
                     
                     count_zeta = count_zeta + len_monom_all;
                 end
                 

                 %ship off dual variables for processing in subsystem 
                 obj.sys{i} = obj.sys{i}.dual_process(obj.dual.v, zeta);       
                 
                 
                 %figure out nonnegativity requirement
             end
            
             %nonnegativity of location (not subsystems)
            %initial measure
            if isempty(obj.init)
                nn_init = 0;
            else
                nn_init = obj.dual.gamma - obj.dual.v;
            end
            
            %terminal measure
            if isempty(obj.term)
                nn_term = 0;
            else
                nn_term = obj.dual.v - obj.dual.beta'*obj.objective;
            end
            obj.dual.nn = [nn_init; nn_term];
             
        end        
        
        
        
        %% Sampling
        
        
        function cb = cost_beta(obj, t, x)
            if isempty(obj.objective)
                cb = zeros(size(t));
            else
                cb = eval(obj.dual.beta'*obj.objective, obj.get_vars(), [t; x]);
            end
        end

%         function nn_out = nonneg(obj, t, x, th, w)
%             %nonnegative functions at this location
%             
%             %initial measure
%             if isempty(obj.init)
%                 nn_init = 0;
%             else
%                 nn_init = obj.dual.gamma - obj.dual.v;
%             end
%             
%             %terminal measure
%             if isempty(obj.term)
%                 nn_term = 0;
%             else
%                 nn_term = obj.dual.v - obj.dual.beta'*obj.objective;
%             end
%             
%             %occupation measure
%             %TODO: replace with subsystem call
%             
% %             nn_occ = -obj.dual.Lv;
% %             if ~isempty(obj.dual.zeta)
% %                 %handle the boxes
% %                 nn_occ = nn_occ + sum(obj.dual.zeta);
% %             end
% %             
% 
%         
%             
%             %box occupation measure
%             
%             
%             %box complement
%             
%             
%             nn = [nn_init; nn_term];
%             
%             
%             % TODO: set to private nn for fast evaluation
%             if obj.loc_supp.TIME_INDEP
%                 nn_out = eval(nn, [obj.x; obj.th; obj.vars.w], [x; th; w]);                                   
%             else
%                 nn_out = eval(nn, obj.get_vars(),  [t; x; th; w]);                                   
%             end
%         end
        
        
        %something about processing dual_rec to get nonnegative functions
        
        

        
%         %% Sampler
%         
%         %TODO: completely rework this section. Use no-class sampler code as
%         %a model for continuous and discrete sampling
%         function out_sim = sample_traj_loc(obj, t0, x0, Tmax, curr_event)
%             %SAMPLE_TRAJ_LOC Sample a single trajectory starting at (t0, x0) in
%             %this location. Stop when the the trajectory hits a guard or
%             %strays outside the location's support region
%             %
%             %curr_event handles the event detection for leaving the support
%             %region, and guards if enabled.
%             %
%             %OUTPUT:
%             %out_sim is a struct holding the simulation output: time,
%             %state, objective, and nonnegative functions from the dual
%             %solution of SDP.
%             
%             if nargin < 5
%                 curr_event = @obj.supp_event;
%             end
%             
%             
%             %simulate the trajectory
%             curr_ode_options = odeset('Events',curr_event, 'RelTol', 1e-7, ...
%                                       'AbsTol', 1e-8, 'MaxStep', 0.01);
%         
%             out_sim = struct;
%             [out_sim.t, out_sim.x] = ode15s(@obj.f_eval, [t0, Tmax], x0, curr_ode_options);
%             
%             %evaluate nonnegative functions
%             if obj.dual.solved
%                 out_sim.nonneg = obj.nonneg(out_sim.t', out_sim.x')';
%             end
%             
%             out_sim.objective = obj.obj_eval(out_sim.t', out_sim.x')';
%             out_sim.id = obj.id;
%         end
%                 
    end
end

