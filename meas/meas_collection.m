classdef meas_collection < handle
    %MEAS_COLLECTION An interface for a collection of measures
    %a class that holds a set of measures in a cell array 'meas'
    
    properties
        vars;       %variables in measures (@mpol)
        meas;       %cell array of measures (@meas)
        
        meas_type = @meas_base; %a handle to the type of measure used 
        %by this collection
    end
    
    methods
        function obj = meas_collection(loc_supp, varnames)
            %MEAS_COLLECTION Construct an instance of this class
            %   Detailed explanation goes here

            %copy over variables specified in varnames
            
            if nargin < 2
                varnames = fields(loc_supp.vars);
            end

            for i = 1:length(varnames)
                curr_var = varnames{i};
                obj.vars.(curr_var) = loc_supp.vars.(curr_var);
            end
        end       
        
        %% support getters
        function supp_out = supp(obj)
            %get the support of all measures
            supp_out = [];
            for i = 1:length(obj.meas)
                supp_out = [supp_out; obj.meas{i}.supp];
            end
        end                
        
        
        %% moment getters
        
        function mass_out = mass(obj)
            %MASS return the mass (moment of 1) of the measure           
            mass_out = 0;
            for i = 1:length(obj.meas)
                mass_out = mass_out + obj.meas{i}.mass();
            end
        end 
        
        function mmmon_out = mom_monom(obj, dmin, dmax)
            %MOM_MONOM moments of monomials
            if nargin < 3
                dmax = dmin;
                dmin = 0;
            end
            
            mmmon_out = 0;
            for i = 1:length(obj.meas)
                mmmon_out = mmmon_out + obj.meas{i}.mom_monom(dmin, dmax);
            end            
        end     
        
        
        function f_new = var_sub_mom(obj, vars_old, f_old)
            %VAR_SUB_MOM returns the moment of f_old with respect to all
            %measures in this terminal set union
            f_new = 0;
            for i = 1:length(obj.meas)
                f_new = f_new + mom(obj.meas{i}.var_sub_end(vars_old, f_old));
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
            %only moments of order 0-2 in square matrix
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
        
        %% variable definition
        %used in meas_def to define a set of new variables
        
        function meas_new = meas_def(obj, varnames, suffix, supp_ref)           
            %MEAS_DEF Define the measures in the collection
            %declare a variable for each measure (index ind in the union)         

            vars_new = struct;
            old_stack =[];
            new_stack = [];

            for i = 1:length(varnames)
                curr_name = varnames{i};                
                
                if isfield(obj.vars, curr_name) 
                    if isempty(obj.vars.(curr_name))
                        vars_new.(curr_name) = [];
                    else
                        curr_var = obj.vars.(curr_name);
                        %declare a new variable
                        new_name = [curr_name, suffix];
                        mpol(new_name, length(curr_var), 1);
                        %load the new variable into vars_new
                        vars_new.(curr_name) = eval(new_name);

                        old_stack = [old_stack; reshape(obj.vars.(curr_name), [], 1)];
                        new_stack = [new_stack; reshape(vars_new.(curr_name), [], 1)];
                    end
                end
            end
                       
            %get the support
            supp_new = subs_vars(supp_ref, old_stack, new_stack);
           
            
            %define the measure
            meas_new = obj.meas_type(vars_new, supp_new);
        end        
    end
end

