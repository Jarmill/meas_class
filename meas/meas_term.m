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
                obj.meas{i} = obj.meas_def({'t', 'x', 'th'}, suffix, supp_curr);
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
                   
        
        %% moment recovery        
        function d_out = mmat(obj)
            %return moment matrix evaluated at current solution
            d_out = cellfun(@(m) m.mmat(), obj.meas, 'UniformOutput', 'false');
%             d_out = double(mmat(obj.meas));
        end
        
        function d_out = mmat_corner(obj)
            %return top-corner moment matrix evaluated at current solution
            %only moments of order 0-2
%             monom_curr = obj.monom(0, 1);
%             mmat_curr = mom(monom_curr*monom_curr');            
%             d_out = double(mmat_curr);
            d_out = cellfun(@(m) m.mmat_corner(), obj.meas, 'UniformOutput', 'false');
        end
        
                
        %% overloads
        function e = isempty(obj)
            %is the support empty?
            %as in supp = []. The harder question would be 'does the basic
            %semialgebraic set formed by the constraints satisfy a
            %nullstellensatz?'
            e = isempty(obj.meas);
        end
        
        
        
        
    end
end

