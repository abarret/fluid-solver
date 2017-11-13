function u = fmg_wrapper(b)
tol = 1.0e-5;
nMax = 5;
global dx; global dy;
b = reshape(b,sqrt(length(b)),sqrt(length(b)));
u = 0*b;
u = fmg(u,b,tol,nMax,dx,dy);
u = u(:);
end