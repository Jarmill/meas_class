classdef subsystem
    %SUBSYSTEM A subsystem x'=f(t, x, th, w, b) of a possibly uncertain 
    %dynamical system in measure analysis    
    
    properties
        vars = struct('t', [], 'x', [], 'th', [], 'w', [], 'b', []);
        
        %measures
        meas_occ = [];  %occupation measure
        meas_box = {};  %box occupation measures
        meas_comp = {}; %box-complement occupation measures                
        
        f = [];
        X = [];
        param = [];
        disturb = [];
        
%         id_str = '';
        id = [];
    end
    
    methods
        function obj = subsystem(loc_id, sys_id, loc_supp)
            %SUBSYSTEM Construct an instance of this class
            
            %TODO: Finish this
            obj.id = [loc_id, sys_id];
        end
        
        
        
        
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

