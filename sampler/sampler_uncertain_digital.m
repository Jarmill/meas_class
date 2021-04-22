classdef sampler_uncertain_digital < sampler_uncertain_interface
    %SAMPLER_UNCERTAIN_DIGITAL Sample in location with uncertainty in
    %discrete-time dynamics
    
   
    
    methods
        function obj = sampler_uncertain_digital(loc,sampler)
            %LOC_SAMPLER Construct an instance of this class
            %   Detailed explanation goes here

%             obj@sampler_interface
            obj@sampler_uncertain_interface(loc, sampler);
            

        end
       
        
        function out_sim = sample_traj(obj,t0, x0, th0, Tmax)
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
            Nt = Tmax - t0;
%             T = (t0:Tmax);
            X = zeros(Nt+1, length(x0));
            W = zeros(Nt+1, length(obj.loc.vars.w));
            X(1, :) = x0;
            x_curr = x0;                 
            
            system_choice = zeros(Nt, 1);
            
            %main solving loop
            k = 0;
            nsys = length(obj.loc.sys);
            while k < Nt
%             for k = 1:Nt
                %find which systems are valid
                possible_sys = [];
                for i = 1:nsys        
                    event_value = obj.loc.sys{i}.supp_eval(k-1, x_curr);

                    if event_value == 1
                        possible_sys = [possible_sys; i];
                    end
                end
                N_possible = length(possible_sys);
                                
                if N_possible == 0
                    break
                end
                
                curr_sys_ind = randi([1, N_possible], 1, 1);
    
                %choose which system to use. curr_event not needed, validity evaluated
                %on next state iteration
                curr_sys = obj.loc.sys{curr_sys_ind};          
                system_choice(k+1) = curr_sys_ind;
                
                %sample the time-varying disturbance
                w_curr = obj.sampler.w();
                %TODO: ensure that 'time' plays nicely
                curr_f = @(t, x) curr_sys.f_eval([x; th0; w_curr]);
                
                %take the step and store results
                x_next = curr_f(k+t0, x_curr);
                X(k+2, :) = x_next;                
                
                if ~isempty(w_curr)
                    W(k+1, :) = w_curr;
                end
                k = k+1;
                x_curr = x_next;

                
            end

            T = t0+(0:k);
            
            %pad an additional disturbance
            if ~isempty(w_curr)
                W(k+1, :) = w_curr;
            end

            %package up output
            out_sim = struct;

            %trajectories
            out_sim.t = T;
            out_sim.x = X;
            out_sim.th = th0;
            out_sim.w = W;
            
            out_sim.break_sys = system_choice;
            
            %TODO: nonnegative evaluation
            out_sim.objective = obj.loc.obj_eval(out_sim.t', out_sim.x')';                                                                
        end
        
        function nonneg_eval =  nonneg_traj(obj, t, x, th, w)
            %NONNEG_TRAJ evaluate nonnegative trajectories 
            %TODO: handle nonnegative trajectories
            
            %remember to scale time as t/obj.Tmax
            
            %evaluate nonnegative functions elsewhere?
            
            nonneg_eval = 0;
            
            
        end
    end
end

