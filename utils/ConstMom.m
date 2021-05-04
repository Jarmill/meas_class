function moments = ConstMom(d, t_span,x0, YalmipBasis)
%CONSTMOM moments of the graph of the function 
%   (t, f(t)) = (t, x0) between the times in tspan
%   up to order d
%   product of lebesgue_t times delta_x

    if length(t_span) == 1;
        t_span = [0; t_span];
    end

    if(~exist('x0','var') || isempty(x0))
        x0 = [];
    end

    %copied from getLebesgueMomentsNew.m (LebesgueBoxMom.m)        
    if(~exist('YalmipBasis','var') || isempty(YalmipBasis))
        YalmipBasis = 0;
    end
    n = 1 + size(x0,1);
 
    if(YalmipBasis == 1)
%         disp('Generating moments in Yalmip basis')
        dv = monpowers(n,d);
    else
%         disp('Generating moments in Gloptipoly basis')
        dv = genPowGlopti(n,d);
    end
%     moments = zeros(size(dv,1),1);
    
    %end copy
    
    %lebesgue distribution in t
    t_power = dv(:, 1);
    
    t_exp = (t_span(2).^(t_power+1) - t_span(1).^(t_power+1))./(t_power+1);
    
    %dirac delta in x    
    if isempty(x0)
        x_exp = 1;
    else
        x_exp = prod(x0'.^dv(:, 2:end), 2);
    end
    
    moments = t_exp .* x_exp;
    
    
    
%     for i = 1:numel(moments)
%         moments(i) = prod((box(2,:).^(dv(i,:)+1) - box(1,:).^(dv(i,:)+1)) ./ (dv(i,:)+1));
%     end
end

