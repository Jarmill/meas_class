classdef meas_term < handle
    %MEAS_TERM A container of terminal measures
    %   realizes unions of terminal measures for a multi-part X_term 
    
    properties
       
        %variables
        vars = struct('t', [], 'x', [], 'th', []);
        
        %measure of variables
        meas = [];
        NT;
        
        %support of measures
        XT = [];
        THT = [];      
    end
    
    methods
        
        %% constructor
        function obj = meas_term(loc_supp, loc_id)
            %MEAS_INIT Construct a measure
            %include the variables and the support         
            
            if nargin < 2
                loc_id = [];
            end
            
            %copy over variables 
%             obj.vars.t  = loc_supp.vars.t;
%             obj.vars.x  = loc_supp.vars.x;
%             obj.vars.th = loc_supp.vars.th;
            varnames = fields(loc_supp.vars);
            for i = 1:length(varnames)
                curr_var = varnames{i};
                obj.vars.(curr_var) = loc_supp.vars.(curr_var);
            end
            
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
                obj.meas{i} = obj.meas_def(suffix, supp_curr);
            end
            
            obj.NT = NT;
            obj.XT = XT_cell;
            obj.THT = THT_cell;                                 
        end

        
        
        
        %monomials: gloptipoly cannot add together monomials from different
        %measures
                
        
        %% measures
        function meas_new = meas_def(obj, suffix, supp_ref)           
            %declare a variable for each measure (index ind in the union)
            vars_new = struct('t', [], 'x', [], 'th', []);           
            varnames = fields(vars_new);
            for i = 1:length(varnames)
                curr_name = varnames{i};
                curr_var = obj.vars.(curr_name);
                
                if ~isempty(curr_var)
                    %declare a new variable
                    new_name = [curr_name, suffix];
                    mpol(new_name, length(curr_var), 1);
                    %load the new variable into vars_new
                    vars_new.(curr_name) = eval(new_name);
                end
%                 obj.vars.(curr_var) = vars.(curr_var);
            end
            
           
                %create new support as well
%                 supp_ref = ;
                supp_new = subs_vars(supp_ref, [obj.vars.t; obj.vars.x; obj.vars.th], ...
                                [vars_new.t; vars_new.x; vars_new.th]);
           
            
            %define the measure
            meas_new = meas_uncertain(vars_new, supp_new);
        end
        
        %% support
        function supp_out = supp(obj)
            %get the support of all measures
            supp_out = [];
            for i = 1:obj.NT
                supp_out = [supp_out; obj.meas{i}.supp];
            end
        end
        
        %% moments
        
        function mass_out = mass(obj)
            %MASS returns the mass (moment of 1) of the measure           
            mass_out = 0;
            for i = 1:obj.NT
                mass_out = mass_out + obj.meas{i}.mass();
            end
        end   
        
        function f_new = var_sub_mom(obj, vars_old, f_old)
            %VAR_SUB_MOM returns the moment of f_old with respect to all
            %measures in this terminal set union
            f_new = 0;
            for i = 1:obj.NT
                f_new = f_new + mom(obj.meas{i}.var_sub_end(vars_old, f_old));
            end
        end   
        
        function mmmon_out = mom_monom(obj, dmin, dmax)
            %MOM_MMON moments of monomials
            if nargin < 3
                dmax = dmin;
                dmin = 0;
            end
            
            mmmon_out = 0;
            for i = 1:obj.NT
                mmmon_out = mmmon_out + obj.meas{i}.mom_monom(dmin, dmax);
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

