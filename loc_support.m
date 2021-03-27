classdef loc_support
    %LOC_SUPPORT Support of a location of a system
    %   Detailed explanation goes here
    
    properties
        %variables of this location
        vars = struct('t', [], 'x', [], 'theta', [], 'w', [], 'b', []);
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
        function obj = loc_support(vars)
            %LOC_SUPPORT Construct an instance of this class
            %   Detailed explanation goes here
            obj.vars = vars;
        end
        
        
        %get time sets

                
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
            if isempty(obj.X_term)
                X_sys = obj.X;
            else
                X_sys = X_sys_in;
            end
        end
        
        function supp_sys = supp_sys_pack(obj, X_sys_in)
            supp_sys = [obj.get_t_supp_sys();
                    X_sys_in;
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

        
        
        
        
%         function supp_out = get_supp(obj, include_disturb)
%             if nargin == 1
%                 include_disturb = 0;
%             end
%             
%             if obj.TIME_INDEP
%                 t_supp = 
%         end
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

