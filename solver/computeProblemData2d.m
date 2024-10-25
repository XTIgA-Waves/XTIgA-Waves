function problemData = computeProblemData2d(sym_u_ex,sym_c,theta)

syms x y t;

c=matlabFunction(sym_c,'Vars',[x,y,t]);
problemData.c=@(x,y,t) c(x,y,t)+0*x+0*y+0*t;

u_ex=matlabFunction(sym_u_ex,'Vars',[x,y,t]);
problemData.u_ex=@(x,y,t) u_ex(x,y,t)+0*x+0*y+0*t;

sym_dx_u_ex=diff(sym_u_ex,x);
dx_u_ex=matlabFunction(sym_dx_u_ex,'Vars',[x,y,t]);
dx_u_ex=@(x,y,t) dx_u_ex(x,y,t)+0*x+0*y+0*t;
sym_dy_u_ex=diff(sym_u_ex,y);
dy_u_ex=matlabFunction(sym_dy_u_ex,'Vars',[x,y,t]);
dy_u_ex=@(x,y,t) dy_u_ex(x,y,t)+0*x+0*y+0*t;
sym_dt_u_ex=diff(sym_u_ex,t);
dt_u_ex=matlabFunction(sym_dt_u_ex,'Vars',[x,y,t]);
dt_u_ex=@(x,y,t) dt_u_ex(x,y,t)+0*x+0*y+0*t;
problemData.dt_u_ex=dt_u_ex;

problemData.grad_u_ex=@(x,y,t) cat (1, ...
    reshape (dx_u_ex(x,y,t), [1, size(x)]), ...
    reshape (dy_u_ex(x,y,t), [1, size(x)]),...
    reshape (dt_u_ex(x,y,t), [1, size(x)]));

sym_dxx_u_ex=diff(sym_dx_u_ex,x);
dxx_u_ex=matlabFunction(sym_dxx_u_ex,'Vars',[x,y,t]);
dxx_u_ex=@(x,y,t) dxx_u_ex(x,y,t)+0*x+0*y+0*t;
sym_dxy_u_ex=diff(sym_dx_u_ex,y);
dxy_u_ex=matlabFunction(sym_dxy_u_ex,'Vars',[x,y,t]);
dxy_u_ex=@(x,y,t) dxy_u_ex(x,y,t)+0*x+0*y+0*t;
sym_dyy_u_ex=diff(sym_dy_u_ex,y);
dyy_u_ex=matlabFunction(sym_dyy_u_ex,'Vars',[x,y,t]);
dyy_u_ex=@(x,y,t) dyy_u_ex(x,y,t)+0*x+0*y+0*t;

% We compute the Hessian only for space variables
problemData.hess_u_ex = @(x,y,t) cat (1, ...        
    reshape (dxx_u_ex(x,y,t), [1, size(x)]), ...
    reshape (dxy_u_ex(x,y,t), [1, size(x)]), ...
    reshape (dxy_u_ex(x,y,t), [1, size(x)]), ...
    reshape (dyy_u_ex(x,y,t), [1, size(x)]));

sym_dx_c=diff(sym_c,x);
sym_dy_c=diff(sym_c,y);
sym_f=-sym_c^2*(sym_dxx_u_ex+sym_dyy_u_ex)+diff(diff(sym_u_ex,t),t)-2*sym_c*(sym_dx_c*sym_dx_u_ex+sym_dy_c*sym_dy_u_ex);
f=matlabFunction(sym_f,'Vars',[x,y,t]);
problemData.f=@(x,y,t) f(x,y,t)+0*x+0*t+0*y;

problemData.grad_u_ex_c2=@(x,y,t) cat (1, ...
    reshape (dx_u_ex(x,y,t).*problemData.c(x,y,t).^2, [1, size(x)]), ...
    reshape (dy_u_ex(x,y,t).*problemData.c(x,y,t).^2, [1, size(x)]),...
    reshape (dt_u_ex(x,y,t).*problemData.c(x,y,t).^2, [1, size(x)]));

problemData.dt_u_ex=dt_u_ex;

problemData.c_dt_u_ex=@(x,y,t) dt_u_ex(x,y,t)*c(x,y,t);

problemData.theta=theta;

problemData.gInitDer = dt_u_ex;

sym_dx_dt_u_ex=diff(sym_dt_u_ex,x);
dx_dt_u_ex=matlabFunction(sym_dx_dt_u_ex,'Vars',[x,y,t]);
dx_u_ex=@(x,y,t) dx_dt_u_ex(x,y,t)+0*x+0*y+0*t;
sym_dy_dt_u_ex=diff(sym_dt_u_ex,y);
dy_dt_u_ex=matlabFunction(sym_dy_dt_u_ex,'Vars',[x,y,t]);
dy_dt_u_ex=@(x,y,t) dy_dt_u_ex(x,y,t)+0*x+0*y+0*t;
sym_dt_dt_u_ex=diff(sym_dt_u_ex,t);
dt_dt_u_ex=matlabFunction(sym_dt_dt_u_ex,'Vars',[x,y,t]);
dt_dt_u_ex=@(x,y,t) dt_dt_u_ex(x,y,t)+0*x+0*y+0*t;

problemData.grad_dt_u_ex=@(x,y,t) cat (1, ...
    reshape (dx_dt_u_ex(x,y,t), [1, size(x)]), ...
    reshape (dy_dt_u_ex(x,y,t), [1, size(x)]),...
    reshape (dt_dt_u_ex(x,y,t), [1, size(x)]));

end