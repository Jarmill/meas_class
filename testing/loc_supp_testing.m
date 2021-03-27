mset clear

mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('th', 1, 1)
mpol('w', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;
vars.th = th;
vars.w = w;

X = x.^2 <= 1;
TH = th.^2 <= 1;
W = sum(w.^2) <= 1;

X_sys = {[x.^2 <= 1], [(x-0.5).^2 <= 0.25]};

lsupp = loc_support(vars);
lsupp.X = X;
lsupp.X_sys = X_sys;
lsupp.param = TH;
lsupp.disturb = W;

X0 = lsupp.supp_init();
Xp = lsupp.supp_term();
Xs = lsupp.supp_sys();