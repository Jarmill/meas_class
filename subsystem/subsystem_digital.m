classdef subsystem_digital < subsystem_interface
    %SUBSYSTEM_DIGITAL A subsystem x+=f(t, x, th, w) of a possibly 
    %uncertain discrete dynamical system in measure analysis    
    
    %TODO: Debug and determine if this class is correct
    
    properties
        %maybe create a subsystem_uncertain?
       varnames = {'t', 'x', 'th', 'w'}; 
       
       %TODO: figure out how to incorporate TIME_INDEP into digital
    end
    
    methods
        %% Constructor
        function obj = subsystem_digital(loc_supp, f, sys_id, loc_id)
            %SUBSYSTEM_DIGITAL Construct a digital (possibly uncertain) 
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
                                     
        end
               
        
        %% Constraints        
        function Ay = cons_liou(obj, d)
            %CONS_LIOU Liouville Equation includes a comparision with 
            %forward operator by a pushforward (discrete systems only)
            
            Ay = obj.meas_occ.mom_push(d, obj.vars, obj.f) - ...
                              obj.meas_occ.mom_monom_proj(d);
        end
                      
        
        
        %% Dual Recovery  
        
        function obj = dual_process(obj, v, zeta)
            %DUAL_PROCESS store dual functions and compute nonnegative
            %functions for this digital subsystem
            pushforward = eval(v, obj.vars.x, obj.f);
            Lv = pushforward - v;
            obj.nn_ = -Lv;
        end
        
        %% Evaluation and Sampling
        %TODO: Implement this
    end
end

