%feas_test
%try and connect together two points by a box-constrained input

%% variables
mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('b', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;
vars.b = b;


%% location support 
lsupp = loc_support(vars);
lsupp.Tmax = 5;
%set has two connected components
% lsupp.X_sys = [x.^2 <= 9; x(1)^2 >= 1]; 
% lsupp.X = [x.^2 <= 9];
% lsupp.X_sys = {[x.^2 <= 3; x(1) >= 1], [x.^2 <= 3; x(1) <= -1]};
lsupp.X_sys = [x.^2 <= 4; x(1)^2 >= 0.5^2]; 
lsupp.X = [x.^2 <= 4];
FEAS = 0;

lsupp.X_init = [2; 2]/2;
if FEAS    
    lsupp.X_term = [2; -2]/2;
else
    lsupp.X_term = [-2; 2]/2;
end

f = 2*b-1;

objective = t; %maximum time?
% objective = 0; % feasibility

%% Solve problem
PM = peak_manager(lsupp, f, objective);

% order = 1;
Nd = 6;
status_rec = zeros(Nd, 1);
time_rec = zeros(Nd, 1);
for order = 1:Nd
    d = 2*order;
    [obj_p, mom_con, supp_con, len_liou, len_abscont] = PM.peak_cons(d);
    sol = PM.peak_solve(obj_p, mom_con,supp_con, 1);
    status_rec(order) = sol.status;
    time_rec(order) = sol.obj_rec;
end

%the minimizing hierarchy returns a rising sequence of lower bounds
%If the path is truly feasible, then Tmin = inf.
%
%the hierarchy could return 0.6 <= 0.7 <= 0.75 <= inf (infeasible)
%as a sequence. Only in the 'inf' return could it be trusted for
%disconnectedness (?)
%a return of 0.6 does not verify feasibility
