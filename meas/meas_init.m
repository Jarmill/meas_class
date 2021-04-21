classdef meas_init < handle
    %MEAS_INIT A container of initial measures
    %   realizes unions of initial measures for a multi-part X_init
    
    properties
       
        %variables
        vars = struct('t', [], 'x', [], 'th', []);
        
        %measure of variables
        meas = [];
        N0;
        
        %support of measures
        X0 = [];
        TH0 = [];      
    end
    
    methods
        
        %% constructor
        function obj = meas_init(loc_supp, loc_id)
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
                obj.meas{i} = obj.meas_def(suffix, supp_curr);
            end
            
            obj.N0 = N0;
            obj.X0 = X0_cell;
            obj.TH0 = TH0_cell;
                       
            
%             obj.vars = struct('t', vars.t, 'x', vars.x);        
%             obj.meas = meas([obj.vars.t; obj.vars.x]);            
%             obj.supp = supp;
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
        
        function supp_out = supp(obj)
            %get the support of all measures
            supp_out = [];
            for i = 1:obj.N0
                supp_out = [supp_out; obj.meas{i}.supp];
            end
        end
        
        function mass_out = mass(obj)
            %MASS return the mass (moment of 1) of the measure           
            mass_out = 0;
            for i = 1:obj.N0
                mass_out = mass_out + obj.meas{i}.mass();
            end
        end                         
        
        function mmmon_out = mom_monom(obj, dmin, dmax)
            %MOM_MMON moments of monomials
            if nargin < 3
                dmax = dmin;
                dmin = 0;
            end
            
            mmmon_out = 0;
            for i = 1:obj.N0
                mmmon_out = mmmon_out + obj.meas{i}.mom_monom(dmin, dmax);
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

