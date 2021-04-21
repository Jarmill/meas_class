classdef sampler_base_discrete  < sampler_base_interface
    %SAMPLER_BASE_DISCRETE Sample in a single location with discrete-time 
    %dynamics and no uncertainty        
    
    
    methods
        function obj = sampler_base_discrete(loc,sampler)
            %LOC_SAMPLER Construct an instance of this class
            %   Detailed explanation goes here

            obj@sampler_base_interface(loc, sampler);
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

