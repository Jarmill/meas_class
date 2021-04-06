mset clear
%% variables

mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('th', 1, 1)

vars = struct;
vars.t = t;
vars.x = x;
vars.th = th;
%initial set
C0 = [1.5; 0];
R0 = 0.4;

% X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);
X0 = [1 1 2 0; 1 -1 0 0];
%% location support 

lsupp = loc_support(vars);
lsupp = lsupp.set_box(4);
lsupp.X_init = X0;
lsupp.X_term = {[x'*x <= 1; x'*x >= 0.5]; [x'*x <= 2; x>=1]};
lsupp.Tmax = 10;
lsupp.param = th^2<=1;

MI = meas_init(lsupp);

MI.mom_monom(2)

MT = meas_term(lsupp);
MT.mom_monom(2)