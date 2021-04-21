classdef loc_sampler < handle
    %LOC_SAMPLER Sample in a location
    %   Detailed explanation goes here
    
    properties
        vars;
        
        %samplers (default, return empty array)
        sampler = struct('x',  @() [], ...
                        'th', @() [], ...
                        'w',  @() []);
        
        loc;
                    
        %function handle
        odefcn = @ode15s;
        
        %expected time to switch systems (continuous time)
        %exponential distribution with parameter mu
        mu = 1;
        
        Tmax = 1;
        
        FINE = 1; %low tolerance in ode solver
    end
    
    methods
        function obj = loc_sampler(loc,sampler)
            %LOC_SAMPLER Construct an instance of this class
            %   Detailed explanation goes here

            %bind this sampler to the location
            obj.loc = loc;
            loc.sampler = obj;
            
%             obj.sampler = sampler;
            %iterate over arguments in point sampler functions
            sampler_names = fields(sampler);
            for i = 1:length(sampler_names )
                curr_smp = sampler_names{i};
                obj.sampler.(curr_smp) = sampler.(curr_smp);
            end
            
        end
        
        function [out_sim_multi] = sample_traj_multi(obj, N, Tmax)
            if nargin < 3
                Tmax = obj.loc.Tmax;
            end
            
            if isnumeric(obj.sampler.x)
                %given sample points
                N = size(obj.sampler.x, 2);
                out_sim_multi = cell(N, 1);               
                %parallel code requires splitting off separate objects
                for i = 1:N                    
                    x0  = obj.sampler.x();                    
                    th0 = obj.sampler.th();
                    out_sim_multi{i} = obj.sample_traj(0, x0, th0, Tmax);
                end
                
            else
                %random sample.               
                out_sim_multi = cell(N, 1);
                for i = 1:N                    
                    x0  = obj.sampler.x();                    
                    th0 = obj.sampler.th();
                    out_sim_multi{i} = obj.sample_traj(0, x0, th0, Tmax);
                end
            end
        end
        
        function out_sim = sample_traj(obj,t0, x0, th0, Tmax)
            %SAMPLE_TRAJ Sample a single trajectory starting at (t0, x0) in
            %this location. Stop when the the trajectory hits a guard or
            %strays outside the location's support region
            %
            %curr_event handles the event detection for leaving the support
            %region, and guards if enabled.
            %
            %OUTPUT:
            %out_sim is a struct holding the simulation output: time,
            %state, objective, and nonnegative functions from the dual
            %solution of SDP.
            if nargin < 5
                curr_event = @obj.loc.supp_event;                
            end
            
            %copy over from peak/sample_cont.m
            %main solving loop
            x0_curr = x0;
            time_accum = [];
            x_accum = [];
            time_index= [];

            w_accum = [];
            b_accum = [];

            time_total = t0;
            
            system_choice = [];
            time_breaks = 0;
            k = 1;
            
            while time_total < Tmax
    
                %choose a possible switching system that is admissible for current time/state               
                possible_sys = [];
                for i = 1:length(obj.loc.sys)        
                    event_value = obj.loc.sys{i}.supp_eval(time_total, x0_curr);
                    if event_value == 1
                        possible_sys = [possible_sys; i];
                    end
                end

                N_possible = length(possible_sys);
                if N_possible == 0
                    break
                end

                %choose random switching system from this possible set
                curr_sys_ind = randi([1, N_possible], 1, 1);

                %for how long should this uncertain system be sampled?
                if (length(obj.loc.sys) == 1) && (isempty(obj.loc.vars.b)) && (isempty(obj.loc.vars.w))
                    %is there any switching? If not, sample until Tmax
                    time_track = Tmax;
                else
                    %is there any switching? If not, sample for a random
                    %interval of time
                    time_track = exprnd(obj.mu, 1, 1);
                end
                %do not exceed Tmax
                time_track_trunc = min(time_track, Tmax - time_total);


                curr_sys = obj.loc.sys{curr_sys_ind};
                w_curr = obj.sampler.w();
                b_curr = rand(length(obj.loc.vars.b), 1);

                curr_f = @(t, x) curr_sys.f_eval([t; x; th0; w_curr; b_curr]);
                %currently w is not supported in event function. change this
                curr_event = @(t, x) curr_sys.supp_event(t + time_total, x);


                
                %simulate the current system
                if obj.FINE
                    curr_ode_options =   odeset('Events',curr_event, 'RelTol', 1e-7, ...
                                                  'AbsTol', 1e-8, 'MaxStep', 0.01);
                else
                    curr_ode_options =   odeset('Events',curr_event);
                end

                [time_curr, x_curr] = obj.odefcn(curr_f, [0, time_track_trunc], x0_curr, curr_ode_options);              

                %save current trajectory
                %check indices/dimensions
                x0_curr = x_curr(end, :)';
                x_accum = [x_accum; x_curr];
                time_accum = [time_accum; time_curr + time_total];

                time_total = time_total + time_curr(end);

                %time dependent uncertainty
                if isempty(w_curr)
                    w_accum = [w_accum; zeros(size(x_curr, 1), 0)];
                else
                    w_accum = [w_accum; ones(size(x_curr, 1), 1)* w_curr'];  %general            
                end
                b_accum = [b_accum; ones(size(x_curr, 1), 1)* b_curr'];  %box
                system_choice = [system_choice; curr_sys];  %system switching
                time_breaks = [time_breaks; time_total];

                time_index= [time_index; k*ones(length(time_curr), 1)];
                k = k + 1;
            end
            
%             %simulate the trajectory
%             curr_ode_options = odeset('Events',curr_event, 'RelTol', 1e-7, ...
%                                       'AbsTol', 1e-8, 'MaxStep', 0.01);
%             
            %package up the output            
            out_sim = struct;

            %trajectories
            out_sim.t = time_accum;
            out_sim.x = x_accum;
            out_sim.th = th0;
            out_sim.w = w_accum;
            out_sim.b = b_accum;
          
%             %evaluate nonnegative functions
%             if obj.loc.dual.solved                
%                 out_sim.nonneg = obj.nonneg_traj(time_accum, x_accum, th0, w_accum);
%             end
            
            out_sim.objective = obj.loc.obj_eval(out_sim.t', out_sim.x')';
            out_sim.id = obj.loc.id;                                                                    
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

