mset clear
%% variables
BOX = 1;

mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('th', 1, 1)
mpol('w', 2, 1)


vars = struct;
vars.t = t;
vars.x = x;
vars.th = th;
vars.w = w;


if BOX
    mpol('b', 2, 1)
    vars.b = b;
end

%% location support 
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

%% copying support
mpol('tq', 1, 1)
mpol('xq', 2, 1)
mpol('thq', 1, 1)
mpol('wq', 2, 1)

varsq = struct;
varsq.t = tq;
varsq.x = xq;
varsq.th = thq;
varsq.w = wq;

qsupp = loc_support(varsq, lsupp);

%% testing subsystem
% f = -x.*w + th + b(1)*x(1)^2 - b(2)*th.*x;

f = -x + 1;
sys = subsystem(lsupp, f, 2);

loc = location(lsupp, f, x(1), 1);