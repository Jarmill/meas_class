classdef meas_uncertain < meas_interface
    %MEAS_UNCERTAIN A generic measure with uncertainty
    %   (time, state, time-independent uncertainty, time-dependent uncertainty)
              
    methods
        
        %% constructor
        function obj = meas_uncertain(vars, supp)
            %MEAS_UNCERTAIN Construct a measure with variables possessing
            %uncertainty

            obj@meas_interface(vars, supp);
            
            if isnumeric(obj.supp)
                supp_new = ([vars.t; vars.x; vars.th] == supp);
                obj.supp = supp_new;
            end
            
            
        end
        
        %% monomials                
        function mmon_out = monom_proj(obj, dmin, dmax)
            %MMON monomials of variables of measure
            %from degree dmin to dmax
            if nargin == 2
                dmax = dmin;
                dmin = 0;
            end
                        
            mmon_out = mmon([obj.vars.t; obj.vars.x; obj.vars.th], dmin, dmax);
            
            if isempty(obj.supp)
                %empty support: zero moments
                %may disable this
                mmon_out = zeros(size(mmon_out));
            end
        end          
               
        function f_new = var_sub_end(obj, vars_old, f_old)
            %substitute variables of measures in for f_old            
            f_new = subs_vars(f_old, [vars_old.t; vars_old.x; vars_old.th], ...
                                [obj.vars.t; obj.vars.x; obj.vars.th]);
        end
              
        %% measures

        
        function mmmon_out = mom_monom_proj(obj, dmin, dmax)
            %MOM_MMON moments of monomials excluding w            
            if nargin < 3
                dmax = dmin;
                dmin = 0;
            end
            
            mmmon_out = mom(obj.monom_proj(dmin, dmax));
        end
        
        function mom_out = mom_lie(obj, d, vars_old, f_old, suppress_time)
            %lie moments
            v = obj.monom_proj(d);
            f_curr = obj.var_sub(vars_old, f_old);
            mom_out = mom(diff(v, obj.vars.x)*f_curr);
            
            if ~isempty(obj.vars.t) || (nargin ==5)
                mom_out = mom(diff(v, obj.vars.t)) + mom_out;
            end
        end
        
        function mom_out = mom_push(obj, d, vars_old, f_old)
            %MOM_PUSH pushforward moments v(f(x)) - v(x)
            v = obj.monom_proj(d);
            f_curr = obj.var_sub([vars_old.t; vars_old.x; vars_old.th; vars_old.w], f_old);
            Rv = subs(v, [obj.vars.t; obj.vars.x; obj.vars.th], ...
                [obj.vars.t; f_curr; obj.vars.th]);
            mom_out = mom(Rv);
        end        
        
        
        %% moment recovery
        function [optimal, mom_out, corner] = recover(obj, tol)
            %RECOVER if top corner of the moment matrix is rank-1, then
            %return approximate optimizer
            
            if nargin < 2
                tol = 5e-4;
            end
            
            corner = obj.mmat_corner();
            
            mass_curr= corner(1, 1);
            if mass_curr < tol*1e-3
                %measure is empty
                %nobody home
                optimal = 1;
                mom_out.t = [];
                mom_out.x = [];
                mom_out.w = [];
                corner = zeros(size(corner));
            else
                rankM = rank(corner, tol);            
                optimal = (rankM == 1);

                mom_out.t = corner(2, 1);
                mom_out.x = corner(2+(1:length(obj.vars.x)), 1);
                mom_out.w = corner((2+length(obj.vars.x)) + (1:length(obj.vars.w)), 1);
            end
        end                
        
        
    end
end

