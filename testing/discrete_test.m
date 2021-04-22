A = [-0.3 0.8;
    -0.75 -0.3];
mpol('x', 2, 1);
vars.x = x;

f = A*x;

C0 = [-1.5; 0];
R0 = 0.4;

%C0 = [0.8; 0];
%R0 = 0.1;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

objective = -x(2);

lsupp = loc_support(vars);
lsupp = lsupp.set_box(4);
lsupp.X_init = X0;
lsupp.Tmax = 20;
lsupp.DIGITAL = 1;
lsupp.TIME_INDEP = 1;

PM = peak_manager(lsupp, f, objective);

order = 3;
d = 2*order;
sol = PM.run(order, 10);

% smp = struct('x', @() circle_sample(1)'*R0 + C0);
% 
% LS = sampler_base_discrete(loc, smp);
% 
% traj = LS.sample_traj(0, C0, 20);