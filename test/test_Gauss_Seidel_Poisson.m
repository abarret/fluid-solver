clear; close all;
addpath('..');
u_exact = @(x,y) exp(-100*((x-0.5).^2+(y-0.5).^2))-pi*erf(5)*(erf(5)+erf(15))/400;
f = @(x,y) 400*exp(-100*((x-0.5).^2+(y-0.5).^2)).*(49+100*(x-1).*x+100*(y-1).*y);
NN = [32 64 128 256];
i = 1;
global dx; global dy;
for N = NN
Nx = N; Ny = N;
Lx = 1; Ly = 1;


x_side = linspace(0,Lx,Nx+1); y_side = linspace(0,Ly,Ny+1);
x = 0.5*(x_side(2:end)+x_side(1:end-1));
y = 0.5*(y_side(2:end)+y_side(1:end-1));
dx = x(2)-x(1); dy = y(2)-y(1);
[xx,yy] = meshgrid(x,y);
ff = f(xx,yy);
u = Gauss_Seidel_Poisson(zeros(Nx,Ny),ff,1.0e-5,100000);
figure(1);
surf(xx,yy,u,'edgecolor','none');
error = abs(u-u_exact(xx,yy));
figure(2); clf;
surf(xx,yy,error,'edgecolor','none');
L1_norm(i) = dx*dy*sum(sum(error));
L2_norm(i) = dx*dy*sum(sum(error.^2));
max_norm(i) = max(max(error));
i = i+1;
pause()
end
figure(3); clf;
loglog(1./NN,L1_norm)
