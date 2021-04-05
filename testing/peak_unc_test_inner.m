mset clear
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('th', 1, 1)
mpol('w', 1, 1)


vars = struct;
vars.t = t;
vars.x = x;
vars.th = th;
vars.w = w;

%initial set
C0 = [-1; -1];
R0 = 0.5;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%% location support 
X = x.^2 <= 1;
wmax = 0.2;
thmax = 0.5;

lsupp = loc_support(vars);
% lsupp.X = x.^2<=9;
lsupp = lsupp.set_box([-2, 2; -2, 4]);
lsupp.X_sys = [];
lsupp.X_init = X0;
lsupp.param =  (th^2 <= thmax^2);
lsupp.disturb = (w^2 <= wmax^2);
lsupp.Tmax = 10;
%% testing peak estimation

%dynamics
% w = 0;
f = [-0.5*x(1) - (0.5 + w)*x(2) + 0.5; -0.5*x(2) + 1 + th];
objective = x(1);

PM = peak_manager(lsupp, f, objective);

%generate constraints
order = 4;
d = 2*order;
[obj_p, mom_con, supp_con, len_liou, len_abscont] = PM.peak_cons(d);
sol = PM.peak_solve(obj_p, mom_con,supp_con);