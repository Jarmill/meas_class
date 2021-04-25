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
%                 suffix = ['_', num2str(loc_id), '_', num2str(i), '_init'];
                supp_curr = [obj.vars.t==0; X0_cell{i}; TH0_cell{i}];
                obj.meas{i} = obj.meas_def({'t', 'x', 'th'}, suffix, supp_curr);
            end
            
%             obj.N0 = N0;
%             obj.X0 = X0_cell;
%             obj.TH0 = TH0_cell;
                       
        end

        
        
        
        %monomials: gloptipoly cannot add together monomials from different
        %measures
                
        %% measures
%         function meas_new = meas_def(obj, suffix, supp_ref)           
%             %MEAS_DEF Define the measures in the collection
%             %declare a variable for each measure (index ind in the union)
%             vars_new = struct('t', [], 'x', [], 'th', []);           
%             varnames = fields(vars_new);
%             for i = 1:length(varnames)
%                 curr_name = varnames{i};
%                 curr_var = obj.vars.(curr_name);
%                 
%                 if ~isempty(curr_var)
%                     %declare a new variable
%                     new_name = [curr_name, suffix];
%                     mpol(new_name, length(curr_var), 1);
%                     %load the new variable into vars_new
%                     vars_new.(curr_name) = eval(new_name);
%                 end
% %                 obj.vars.(curr_var) = vars.(curr_var);
%             end
%                        
%                 supp_new = subs_vars(supp_ref, [obj.vars.t; obj.vars.x; obj.vars.th], ...
%                                 [vars_new.t; vars_new.x; vars_new.th]);
%            
%             
%             %define the measure
%             meas_new = meas_uncertain(vars_new, supp_new);
%         end
        
                        
        
                               
        
        
        
    end
end

