classdef meas_base < handle
    %MEAS_BASE A measure used for 
    %   Detailed explanation goes here
    
    properties
       
        %variables
        vars = struct('t', [], 'x', []);
        
        %measure of variables
        meas = [];
        
        %support of measures
        supp = [];
                
    end
    
    methods
        
        %% constructor
        function obj = meas_base(vars, supp)
            %MEAS_BASE Construct a measure
            %include the variables and the support         

            if isfield(vars, 'supp')
                supp = vars.supp;
            end
            
            if isnumeric(supp)
                supp_new = ([vars.t; vars.x] == supp);
                supp = supp_new;
            end
            
            obj.vars = struct('t', vars.t, 'x', vars.x);        
            obj.meas = meas([obj.vars.t; obj.vars.x]);            
            obj.supp = supp;
        end
        
        function vars_out = get_vars(obj)
            %GET_VARS add more variables as necessary
            vars_out = [obj.vars.t; obj.vars.x];
        end
        
        %% monomials        
        function mmon_out = monom(obj, dmin, dmax)
            %MMON monomials of variables of measure
            %from degree dmin to dmax
            if nargin == 2
                dmax = dmin;
                dmin = 0;
            end
                        
            mmon_out = mmon([obj.vars.t; obj.vars.x], dmin, dmax);
            
            if isempty(obj.supp)
                %empty support: zero moments
                %may disable this
                mmon_out = zeros(size(mmon_out));
            end
        end     
        
        
        function f_new = var_sub(obj, vars_old, f_old)
            %substitute variables of measures in for f_old            
            f_new = subs(f_old, [vars_old.t; vars_old.x], ...
                                [obj.vars.t; obj.vars.x]);
        end
                
        %% measures
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
        
        function mom_out = mom_lie(obj, d, vars_old, f_old)
            %lie moments
            v = obj.monom(d);
            f_curr = obj.var_sub(vars_old, f_old);
            mom_out = mom(diff(v, obj.vars.t) + diff(v, obj.vars.x)*f_curr);
        end
        
        function mom_out = mom_push(obj, d, vars_old, f_old)
            %pushforward moments v(f(x)) - v(x)
            v = obj.monom(d);
            f_curr = obj.var_sub(vars_old, f_old);
            Rv = subs(v, [obj.vars.t; obj.vars.x], [obj.vars.t; f_curr]);
            mom_out = mom(Rv);
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
        
        function [optimal, mom_out] = recover(obj, tol)
            %RECOVER if top corner of the moment matrix is rank-1, then
            %return approximate optimizer
            
            if nargin < 2
                tol = 5e-4;
            end
            
            corner = obj.mmat_corner();
            rankM = rank(corner, 1e-3);
            optimal = (rankM == 1);
    
            mom_out.t = double(mom(obj.vars.t));
            mom_out.x = double(mom(obj.vars.x));
    
        end
        
        %% overloads
        function e = isempty(obj)
            %is the support empty?
            %as in supp = []. The harder question would be 'does the basic
            %semialgebraic set formed by the constraints satisfy a
            %nullstellensatz?'
            e = isempty(obj.supp);
        end
        
        
        
        
    end
end

