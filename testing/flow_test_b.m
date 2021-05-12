mset clear
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('b', 1, 1);
vars = struct;
vars.t = t;
vars.x = x;
vars.b = b;

%initial set
C0 = [1.5; 0];
R0 = 0.4;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%% location support 

lsupp = loc_support(vars);
% lsupp = lsupp.set_box(4);
lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
lsupp.X_init = X0;
lsupp.Tmax = 5;
%% testing peak estimation

%dynamics
% f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

% dmax = 0.5;
% dmax = 0.4;
% dmax = 0.3;
% dmax = 0.2;
dmax = 0.15;
draw = dmax.*(2.*b - 1);
% f = [x(2); -x(1)*(1+ draw) - x(2) + (1/3).* x(1).^3 ];
f0 = [x(2); -(1-dmax)*x(1)-x(2)+(x(1)^3)/3];
f1 = [0; -2*dmax*x(1)];
f = f0 + b*f1;

% objective = -x(2);
objective = -sum(x);
% objective = -x;

PM = peak_manager(lsupp, f, objective);

%generate constraints
order = 4;
% d = 2*order;
sol = PM.run(order, lsupp.Tmax);
sol.obj_rec
% [obj_p, mom_con, supp_con, len_liou, len_abscont] = PM.peak_cons(d);
% sol = PM.peak_solve(obj_p, mom_con,supp_con);