classdef sampler_uncertain_interface < handle
    %SAMPLER_UNCERTAIN_INTERFACE Sampling in a location with uncertainty
    %   Detailed explanation goes here
    
    properties
        loc;     %location
        
        vars;    %variables in locations
        sampler = struct('x', @() [], 'w', @() [], 'th', @() []); %samplers (default, return empty array)
                        
        %function handle
        odefcn = @ode15s;
        
        Tmax = 1;
        
        %low tolerance in ode solver
        FINE = 1;
    end
    
    methods
        function obj = sampler_uncertain_interface(loc,sampler)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            
            %bind this sampler to the given location
            obj.loc = loc;
            loc.sampler = obj;
            
            
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
        
    end
    
    methods(Abstract)
        sample_traj(obj, t0, x0)
        %sample a trajectory starting at a single location
        
        nonneg_eval =  nonneg_traj(obj, t, x)
        %evaluate nonnegative functions along trajectories
    end
end

