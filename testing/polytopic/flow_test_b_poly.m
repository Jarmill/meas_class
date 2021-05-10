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
lsupp.Tmax = 10;

% lsupp.poly = struct('A', [1; -1], 'b', [-1; -1]);

%try out A(input) <= b instead
%with a centrally symmetric constraint set
lsupp.poly = struct('A', [1], 'b', [1]);
%% testing peak estimation

%dynamics
% f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

dmax = 0.5;
% dmax = 0.4;
% dmax = 0.3;
% dmax = 0.2;
% draw = dmax.*(2.*b - 1);
draw = dmax * b;
f = [x(2); -x(1) - x(2) + (1/3).* x(1).^3 + draw];

objective = -x(2);

PM = peak_manager(lsupp, f, objective);

%generate constraints
order = 5;
% d = 2*order;
sol = PM.run(order, lsupp.Tmax);
% [obj_p, mom_con, supp_con, len_liou, len_abscont] = PM.peak_cons(d);
% sol = PM.peak_solve(obj_p, mom_con,supp_con);