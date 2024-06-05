function problemData = computeProblemData1d(sym_u_ex,sym_c,theta)

syms x t;

c=matlabFunction(sym_c,'Vars',[x,t]);
problemData.c=@(x,t) c(x,t)+0*x+0*t;

u_ex=matlabFunction(sym_u_ex,'Vars',[x,t]);
problemData.u_ex=@(x,t) u_ex(x,t)+0*x+0*t;

sym_dx_u_ex=diff(sym_u_ex,x);
dx_u_ex=matlabFunction(sym_dx_u_ex,'Vars',[x,t]);
dx_u_ex=@(x,t) dx_u_ex(x,t)+0*x+0*t;
sym_dt_u_ex=diff(sym_u_ex,t);
dt_u_ex=matlabFunction(sym_dt_u_ex,'Vars',[x,t]);
dt_u_ex=@(x,t) dt_u_ex(x,t)+0*x+0*t;
problemData.dt_u_ex=dt_u_ex;

problemData.grad_u_ex=@(x,t) cat (1, ...
    reshape (dx_u_ex(x,t), [1, size(x)]), ...
    reshape (dt_u_ex(x,t), [1, size(x)]));
sym_dxx_u_ex=diff(sym_dx_u_ex,x);
dxx_u_ex=matlabFunction(sym_dxx_u_ex,'Vars',[x,t]);
dxx_u_ex=@(x,t) dxx_u_ex(x,y,t)+0*x+0*t;


% We compute the Hessian only for the space variable
problemData.hess_u_ex = dxx_u_ex;

sym_dx_c=diff(sym_c,x);
sym_f=-sym_c^2*sym_dxx_u_ex+diff(diff(sym_u_ex,t),t)-2*sym_c*sym_dx_c*sym_dx_u_ex;
f=matlabFunction(sym_f,'Vars',[x,t]);
problemData.f=@(x,t) f(x,t)+0*x+0*t;

problemData.grad_u_ex_c2=@(x,t) cat (1, ...
    reshape (dx_u_ex(x,t).*problemData.c(x,t).^2, [1, size(x)]), ...
    reshape (dt_u_ex(x,t).*problemData.c(x,t).^2, [1, size(x)]));

problemData.dt_u_ex=dt_u_ex;

problemData.c_dt_u_ex=@(x,t) dt_u_ex(x,t)*c(x,t);

problemData.theta=theta;

end