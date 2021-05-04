classdef meas_base < meas_interface
    %MEAS_BASE A generic measure with only (time, state)
    %refer to meas_uncertain for uncertainty
    
    methods        
        %% constructor
        function obj = meas_base(vars, supp)
            %MEAS_BASE Construct a measure
            %include the variables and the support         

            obj@meas_interface(vars, supp);
            
            if isnumeric(obj.supp) && ~isempty(obj.supp)
                supp_new = ([vars.t; vars.x] == supp);
                obj.supp = supp_new;
            end
            
            
        end
              
        %% liouville moments 
        function mom_out = mom_lie(obj, d, vars_old, f_old, suppress_time)
            %lie moments
            v = obj.monom(d);
            f_curr = obj.var_sub(vars_old, f_old);
            mom_out = mom(diff(v, obj.vars.x)*f_curr);
            
            if (isfield(obj.vars, 't') && ~isempty(obj.vars.t)) && (nargin < 5)
                mom_out = mom(diff(v, obj.vars.t)) + mom_out;
            end
        end
        
        function mom_out = mom_push(obj, d, vars_old, f_old)
            %MOM_PUSH pushforward moments v(f(x)) - v(x)
            v = obj.monom_proj(d);
            f_curr = obj.var_sub(vars_old, f_old);
            Rv = subs(v, [obj.vars.t; obj.vars.x], ...
                [obj.vars.t; f_curr]);
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
                corner = zeros(size(corner));
            else
                rankM = rank(corner, tol);            
                optimal = (rankM == 1);

                if (isfield(obj.vars, 't') && ~isempty(obj.vars.t))
                    mom_out.t = corner(2, 1);
                    mom_out.x = corner(2+(1:length(obj.vars.x)), 1);                
                else
                    mom_out.x = corner(1+(1:length(obj.vars.x)), 1);                
                end
            end
        end                
        
        
    end
end

