classdef sampler_base < sampler_base_interface
    %SAMPLER_BASE Sample in a single location with continuous-time dynamics 
    %and no uncertainty        
    
    methods
        function obj = sampler_base(loc,sampler)
            %LOC_SAMPLER Construct an instance of this class
            %   Detailed explanation goes here

            obj@sampler_base_interface(loc, sampler);
        end
        
        
        
        function out_sim = sample_traj(obj,t0, x0, Tmax)
            %SAMPLE_TRAJ Sample a single trajectory starting at (t0, x0) in
            %this location. Stop when the the trajectory strays outside the 
            %location's support region
            %
            %curr_event handles the event detection for leaving the support            
            %
            %OUTPUT:
            %out_sim is a struct holding the simulation output: time,
            %state, objective, and nonnegative functions from the dual
            %solution of LMI (TODO).
            
            %only one location
            curr_event = @obj.loc.supp_event;
            
            %dynamics
            curr_sys = obj.loc.sys{1};            
            curr_f = @(t, x) curr_sys.f_eval([t; x]);

            %simulate the current system
            if obj.FINE
                curr_ode_options =   odeset('Events',curr_event, 'RelTol', 1e-7, ...
                                              'AbsTol', 1e-8, 'MaxStep', 0.01);
            else
                curr_ode_options =   odeset('Events',curr_event);
            end

            [time_curr, x_curr] = obj.odefcn(curr_f, [t0, Tmax], x0, curr_ode_options);              
   
            %package up the output            
            out_sim = struct;

            %trajectories
            out_sim.t = time_curr;
            out_sim.x = x_curr;         
            
            out_sim.objective = obj.loc.obj_eval(out_sim.t', out_sim.x')';
            out_sim.id = obj.loc.id;                                                                    
        end
        
        function nonneg_eval =  nonneg_traj(obj, t, x)
            %NONNEG_TRAJ evaluate nonnegative trajectories 
            %TODO: handle nonnegative trajectories
            
            %remember to scale time as t/obj.Tmax
            
            %evaluate nonnegative functions elsewhere?
            
            nonneg_eval = 0;                        
        end
    end
end

