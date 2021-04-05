mset clear
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;

%initial set
C0 = [1.5; 0];
R0 = 0.4;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%% location support 

lsupp = loc_support(vars);
lsupp = lsupp.set_box(4);
lsupp.X_init = X0;
lsupp.Tmax = 10;
%% testing peak estimation

%dynamics
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
objective = -x(2);

PM = peak_manager(lsupp, f, objective);

%generate constraints
order = 3;
d = 2*order;
[obj_p, mom_con, supp_con, len_liou, len_abscont] = PM.peak_cons(d);
sol = PM.peak_solve(obj_p, mom_con,supp_con);