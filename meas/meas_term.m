classdef meas_term < meas_collection
    %MEAS_TERM A container of terminal measures
    %   realizes unions of terminal measures for a multi-part X_term    
    
    methods
        
        %% constructor
        function obj = meas_term(loc_supp, loc_id)
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
            XT = loc_supp.get_X_term();
            
            %NT: number of initial regions
            if isnumeric(XT)
                NT = size(XT, 2);
                XT_cell = cell(NT, 1);
                for i = 1:NT
                    XT_cell{i} = (obj.vars.x == XT(:, i));
                end
                XT = XT_cell;
            else
                if iscell(XT)
                    XT_cell = XT;
                else
                    XT_cell = {XT};
                end
                %XT is a cell
                %could make this less restrictive later
                NT = length(XT_cell);
                
            end
            
            %process disturbance THT
            THT = loc_supp.param;
            if ~isempty(THT)                
                if isnumeric(THT)
                    NTHT = size(THT, 2);
                    THT_cell = cell(NTHT , 1);
                    for i = 1:NTHT
                        THT_cell{i} = (obj.vars.th == THT(:, i));
                    end
                else    
                    if iscell(THT)
                        THT_cell = THT;
                    else
                        THT_cell = {THT};
                    end
                    NTHT = length(THT);
                end
                
                if (NT > 1) && (NTHT == 1)
                    THT_cell2 = cell(NT, 1);
                    for i = 1:NT
                        THT_cell2{i} = THT;
                    end
                    THT_cell = THT_cell2;
                end
            else
                THT_cell = cell(NT, 1);
            end
            
            %now form measures
            obj.meas = cell(NT, 1);
            tsupp =loc_supp.get_t_supp_term();
            for i = 1:NT
                suffix = '_term';
                if NT > 1
                    suffix = ['_', num2str(i), suffix];
                end
                if ~isempty(loc_id)
                    suffix = ['_', num2str(loc_id), suffix];
                end

                supp_curr = [tsupp ; XT_cell{i}; THT_cell{i}];
                obj.meas{i} = obj.meas_def(varnames, suffix, supp_curr);
            end                               
        end                                                          
        
        function mmmon_out = mom_monom_x(obj, dmin, dmax)
            %MMON monomials of variables of measure
            %from degree dmin to dmax
            if nargin == 2
                dmax = dmin;
                dmin = 0;
            end
                
            mmmon_out = 0;
            for i = 1:length(obj.meas)
                mmon_out = mmon(obj.meas{i}.vars.x, dmin, dmax);
                mmmon_out = mmmon_out + mom(mmon_out);
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

