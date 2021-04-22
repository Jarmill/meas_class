classdef sampler_base_discrete  < sampler_base_interface
    %SAMPLER_BASE_DISCRETE Sample in a single location with discrete-time 
    %dynamics and no uncertainty        
    
    
    methods
        function obj = sampler_base_discrete(loc,sampler)
            %LOC_SAMPLER Construct an instance of this class
            %   Detailed explanation goes here

            obj@sampler_base_interface(loc, sampler);
        end
        
        function out_sim = sample_traj(obj,t0, x0, Tmax)
            %SAMPLE_TRAJ Sample a single trajectory starting at (t0, x0, th0) 
            %in this location. Stop when the the trajectory hits strays 
            %outside the location's support region
            %
            %curr_event handles the event detection for leaving the support
            %region, and guards if enabled.
            %
            %OUTPUT:
            %out_sim is a struct holding the simulation output: time,
            %state, objective, and nonnegative functions from the dual
            %solution of SDP.
            
%             outputArg = obj.Property1 + inputArg;

            Nt = Tmax - t0;
%             T = (t0:Tmax);
            X = zeros(Nt+1, length(x0));
            
            X(1, :) = x0;
            x_curr = x0;
            
            %only one location
            curr_event = @(t, x) obj.loc.supp_eval(t, x);
            
            %dynamics
            curr_sys = obj.loc.sys{1};            
            curr_f = @(t, x) curr_sys.f_eval([x]);

            
            %main solving loop
            k = 0;
            while k < Nt
%             for k = 1:Nt
                
                event_value = curr_event(k+t0+1, x_curr);
                if ~event_value
                    break
                end
                
                x_next = curr_f(k+t0, x_curr);
                X(k+2, :) = x_next;                
                k = k+1;
                x_curr = x_next;
            end
            
            T = t0+(0:k);

            %package up output
            out_sim = struct;

            %trajectories
            out_sim.t = T;
            out_sim.x = X;
            
            %TODO: nonnegative evaluation
            out_sim.objective = obj.loc.obj_eval(out_sim.t', out_sim.x')';
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

