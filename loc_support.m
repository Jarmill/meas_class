classdef loc_support
    %LOC_SUPPORT Support of a location of a system
    %   Detailed explanation goes here
    
    properties
        %variables of this location
        vars = struct('t', [], 'x', [], 'th', [], 'w', [], 'b', []);
        %    t:  time
        %    x:  state
        %theta:  time-independent uncertainty
        %    w:  time-dependent uncertainty
        %    b:  time-dependent box uncertainty in [0, 1] (if not digital)
        
        %uncertainty could also be input, like for control
        
        %time
        Tmax(1,1) double{mustBePositive}  = 5;
        
        
        TIME_INDEP = 0; %
        FREE_TERM = 1;  %free terminal time between 0 and Tmax
        DIGITAL = 0;
        
        %state:
        %all space
        X = [];
        
        %state initial support
        X_init = [];
        
        %state terminal support
        X_term = [];
        
        %state occupation support (subsystems if switching is allowed)
        X_sys = [];
        
        % time-dependent uncertainty
        disturb = []; %(w)
        
        % time-independent uncertainty
        param = []; %(theta)
    end
    
    methods
        
        %% constructor
        function obj = loc_support(vars, loc_ref, sys_id)
            %LOC_SUPPORT Construct an instance of this class
            %INPUT: 
            %   vars:   variables to define support
            %   all other attributes are defined outside the constructor
            %
            %   REFERENCE:  substitute in quantities in reference with
            %               variables in vars
            %   loc_ref:    reference location
            %   sys_id:     system in reference to examine
            
            
            %iterate through variables 
            varnames = fields(vars);
            for i = 1:length(varnames)
                curr_var = varnames{i};
                obj.vars.(curr_var) = vars.(curr_var);
            end
            
            if nargin > 1
                %substitue all attributes of reference with new variables
                %a constructor with copying and substitution
                
                %constants
                obj.Tmax = loc_ref.Tmax;
                obj.TIME_INDEP = loc_ref.TIME_INDEP; %
                obj.FREE_TERM =  loc_ref.FREE_TERM;  %free terminal time between 0 and Tmax
                obj.DIGITAL = loc_ref.DIGITAL; %
                
                %support sets
                obj.X = subs_vars(loc_ref.X, loc_ref.vars.x, obj.vars.x);
                obj.X_init = subs_vars(loc_ref.X_init, loc_ref.vars.x, obj.vars.x);
                obj.X_term = subs_vars(loc_ref.X_term, loc_ref.vars.x, obj.vars.x);
                
                if iscell(loc_ref.X_sys) && (nargin == 2)
                    %copy all subsystems
                    obj.X_sys = cell(length(loc_ref.X_sys), 1);
                    for i = 1:length(loc_ref.X_sys)
                        obj.X_sys{i} = subs_vars(loc_ref.X_sys{i}, loc_ref.vars.x, obj.vars.x);
                    end
                else
                    %copy only one subsystem
                    if nargin == 3
                        %only the subsystem required
                        obj.X_sys = subs_vars(loc_ref.X_sys{i}, loc_ref.vars.x, obj.vars.x);
                    else
                        obj.X_sys = subs_vars(loc_ref.X_sys, loc_ref.vars.x, obj.vars.x);
                    end
                end
                
                obj.disturb = subs_vars(loc_ref.disturb, loc_ref.vars.w, obj.vars.w);
                
                obj.param = subs_vars(loc_ref.disturb, loc_ref.vars.th, obj.vars.th);
                
            end
        end
        
        %% get variables
        
        function vars_out = get_vars(obj)
            vars_out = [obj.vars.t;
                        obj.vars.x;
                        obj.vars.th;
                        obj.vars.w];
        end
        
        
        function vars_out = get_vars_box(obj)
            vars_out = [obj.vars.t;
                        obj.vars.x;
                        obj.vars.th;
                        obj.vars.w;
                        obj.vars.b];
        end
        
        %% get initial set
        
        function t_supp = get_t_supp_init(obj)
            if obj.TIME_INDEP
                t_supp =[];
            else
                t_supp = obj.vars.t == 0;
            end
        end 
        
        function X_init = get_X_init(obj)
            if isempty(obj.X_init)
                X_init = obj.X;
            else
                X_init = obj.X_init;
            end
        end
        
        function supp_out = supp_init(obj)
            %initial set 
            supp_out = [obj.get_t_supp_init();
                        obj.get_X_init();
                        obj.param];                    
        end
        
        
        %% get terminal set                 
        function t_supp = get_t_supp_term(obj)
            if obj.TIME_INDEP
                t_supp =[];
            else
                if obj.FREE_TERM %free terminal time
                    t_supp = obj.vars.t*(obj.Tmax - obj.vars.t)>= 0;
                else
                    t_supp = (obj.vars.t == obj.Tmax);
                end
            end
        end 
        
        function X_term = get_X_term(obj)
            if isempty(obj.X_term)
                X_term = obj.X;
            else
                X_term = obj.X_term;
            end
        end
        
        function supp_out = supp_term(obj)
            %initial set 
            supp_out = [obj.get_t_supp_term();
                        obj.get_X_init();
                        obj.param];                    
        end
        
        %% get system set
        
        function t_supp = get_t_supp_sys(obj)
            if obj.TIME_INDEP
                t_supp =[];
            else
                t_supp = obj.vars.t*(obj.Tmax - obj.vars.t)>= 0;
            end
        end 
        
        function X_sys = get_X_sys_single(obj, X_sys_in)
            if isempty(X_sys_in)
                X_sys = obj.X;
            else
                X_sys = X_sys_in;
            end
        end
        
        function supp_sys = supp_sys_pack(obj, X_sys_in)
            %assumptions: systems are valid for all time
            supp_sys = [obj.get_t_supp_sys();
                    obj.get_X_sys_single(X_sys_in);
                    obj.param;
                    obj.disturb];    
        end
        
        %TODO: deal with switching, where X has multiple sets
        %will need to make a class called 'subsystem'
        function supp_out = supp_sys(obj)
            %initial set 
            
            if iscell(obj.X_sys)
%                 Nsys = length(obj.X_sys);
                supp_out = cellfun(@(Xs) obj.supp_sys_pack(Xs), obj.X_sys,...
                    'UniformOutput', false);
%                 supp_out = cell(Nsys, 1);
%                 for i = 1:Nsys
%                    supp_out{i} =  
%                 end                
            else
                supp_out = obj.supp_sys_pack(obj.X_sys);

            end                            
        end
        
        
                              
    end
end

