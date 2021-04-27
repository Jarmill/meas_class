classdef meas_init < meas_collection
    %MEAS_INIT A container of initial measures
    %   realizes unions of initial measures for a multi-part X_init
    
    methods
        
        %% constructor
        function obj = meas_init(loc_supp, loc_id)
            %MEAS_INIT Construct a measure
            %include the variables and the support         
            
            if nargin < 2
                loc_id = [];
            end

            %copy over variables 
            varnames = {'t','x','th'};
            obj@meas_collection(loc_supp, varnames);            
            obj.meas_type = @meas_uncertain;
            
            %process initial region
            X0 = loc_supp.get_X_init();
            
            %N0: number of initial regions
            if isnumeric(X0)
                N0 = size(X0, 2);
                X0_cell = cell(N0, 1);
                for i = 1:N0
                    X0_cell{i} = (obj.vars.x == X0(:, i));
                end
                X0 = X0_cell;
            else
                if iscell(X0)
                    X0_cell = X0;
                else
                    X0_cell = {X0};
                end
                %X0 is a cell
                %could make this less restrictive later
                N0 = length(X0_cell);
                
            end
            
            %process disturbance TH0
            TH0 = loc_supp.param;
            if ~isempty(TH0)                
                if isnumeric(TH0)
                    NTH0 = size(TH0, 2);
                    TH0_cell = cell(NTH0 , 1);
                    for i = 1:NTH0
                        TH0_cell{i} = (obj.vars.th == TH0(:, i));
                    end
                else    
                    if iscell(TH0)
                        TH0_cell = TH0;
                    else
                        TH0_cell = {TH0};
                    end
                    NTH0 = length(TH0);
                end
                
                if (N0 > 1) && (NTH0 == 1)
                    TH0_cell2 = cell(N0, 1);
                    for i = 1:N0
                        TH0_cell2{i} = TH0;
                    end
                    TH0_cell = TH0_cell2;
                end
            else
                TH0_cell = cell(N0, 1);
            end
            
            %now form measures
            obj.meas = cell(N0, 1);
            for i = 1:N0
                suffix = '_init';
                if N0 > 1
                    suffix = ['_', num2str(i), suffix];
                end
                if ~isempty(loc_id)
                    suffix = ['_', num2str(loc_id), suffix];
                end
                supp_curr = [obj.vars.t==0; X0_cell{i}; TH0_cell{i}];
                obj.meas{i} = obj.meas_def({'t', 'x', 'th'}, suffix, supp_curr);
            end
                                   
        end

        %% recovery
        function [optimal, mom_out, corner] = recover(obj, tol)
            %RECOVER if top corner of the moment matrix is rank-1, then
            %return approximate optimizer
            if nargin < 2
                tol = 5e-4;
            end
            
            Nmeas = length(obj.meas);
            optimal_cell = zeros(Nmeas, 1);
            mom_out_cell = cell(Nmeas, 1);
            corner_cell  = cell(Nmeas, 1);
            for i = 1:Nmeas
                [optimal_cell(i), mom_out_cell{i}, corner_cell{i}] = obj.meas{i}.recover(tol);
            end
            
            if any(optimal_cell)
                optimal = 1;
                mom_out = mom_out_cell{optimal_cell == 1};
                corner= corner_cell{optimal_cell == 1};
            else
                optimal = 0;
                mom_out = struct('t', [], 'x', [], 'th', []);
                corner = [];
            end
        end
        
    end
end

