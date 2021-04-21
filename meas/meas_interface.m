classdef meas_interface < handle
    %MEAS_INTERFACE An interface (abstract) method containing behaviours of
    %measures. These measures can be genericized (such as initial, peak,
    %terminal, occupation measures, etc.)        
    
    properties
        vars;
        meas;
        supp;
    end
    
    methods
        function obj = meas_interface(vars,supp)
            %MEAS_BASE Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;

            if isfield(vars, 'supp')
                obj.supp = vars.supp;
            else
                obj.supp = supp;
            end
            
            varnames = fields(vars);
            for i = 1:length(varnames)
                curr_var = varnames{i};
                obj.vars.(curr_var) = vars.(curr_var);
            end
            
            obj.meas = meas(obj.get_vars());            
            
        end
        
        %% Getters
        function vars_out = get_vars(obj)
            %get variables in measure
            varnames = fields(obj.vars);
            vars_out = [];
            for i = 1:length(varnames)
                curr_var = varnames{i};
                vars_out = [vars_out; obj.vars.(curr_var)];
            end
        end
        
        %% Monomials and Polynomials
        function mmon_out = monom(obj, dmin, dmax)
            %MMON monomials of variables of measure
            %from degree dmin to dmax
            if nargin == 2
                dmax = dmin;
                dmin = 0;
            end
                        
            mmon_out = mmon(obj.get_vars(), dmin, dmax);
            
            if isempty(obj.supp)
                %empty support: zero moments
                %may disable this
                mmon_out = zeros(size(mmon_out));
            end
        end    
            
        function f_new = var_sub(obj, vars_old, f_old)
            %substitute variables of measures in for f_old   
                varnames = fields(vars_old);
                old_stack = [];
                for i = 1:length(varnames)
                    curr_var = varnames{i};
                    old_stack = [old_stack; vars_old.(curr_var)];
                end
            
                f_new = subs_vars(f_old, old_stack, obj.get_vars);
        end
        
        %% Measures and Moments
        function mass_out = mass(obj)
            %MASS return the mass (moment of 1) of the measure           
            if isempty(obj.supp)
                mass_out = 0;
            else
                mass_out = mass(obj.meas);
            end
        end                         
        
        function mmmon_out = mom_monom(obj, dmin, dmax)
            %MOM_MMON moments of monomials
            if nargin < 3
                dmax = dmin;
                dmin = 0;
            end
            
            mmmon_out = mom(obj.monom(dmin, dmax));
        end
        
        function e = isempty(obj)
            %is the support empty?
            %as in supp = []. The harder question would be 'does the basic
            %semialgebraic set formed by the constraints satisfy a
            %nullstellensatz?'
            e = isempty(obj.supp);
        end
        
        %% moment recovery
        
        function d_out = mmat(obj)
            %return moment matrix evaluated at current solution
            d_out = double(mmat(obj.meas));
        end
        
        function d_out = mmat_corner(obj)
            %return top-corner moment matrix evaluated at current solution
            %only moments of order 0-2
            monom_curr = obj.monom(0, 1);
            mmat_curr = mom(monom_curr*monom_curr');            
            d_out = double(mmat_curr);
        end
        
        function d_out = mvec(obj)
            %return moment sequence evaluated at current solution
            d_out = double(mvec(obj.meas));
        end
        
    end
    
    %% Overloads by inheritence
    methods(Abstract)
        
        mom_lie(obj, d, vars_old, f_old, suppress_time)
        %moments of lie derivative
        
        mom_push(obj, d, vars_old, f_old)
        %pushforward moments v(f(x)) - v(x)
                
        recover(obj, tol)
        %if top corner of the moment matrix is rank-1, then
        %return approximate optimizer
    end
end

