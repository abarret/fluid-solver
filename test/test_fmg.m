clear; close all;
addpath('..');
%u_exact = @(x,y) exp(-100*((x-0.5).^2+(y-0.5).^2))-pi*erf(5)*(erf(5)+erf(15))/400;
%f = @(x,y) 400*exp(-100*((x-0.5).^2+(y-0.5).^2)).*(49+100*(x-1).*x+100*(y-1).*y);
u_exact = @(x,y) sin(2*pi*x).*sin(2*pi*y);
f = @(x,y) -8*pi*pi*u_exact(x,y);
NN = [64 128 256 512 1024 2048];
tol = 1.0e-5; 
Nmax = 3;
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
    i = 1;
    [u,b] = fmg(zeros(Nx,Ny),ff,tol,Nmax,dx,dy);
    err(i) = dx*dy*sqrt(sum(sum(b.^2)));
    while( err(i) > tol )
        i = i+1;
        [u,b] = fmg(u,ff,tol/100,Nmax,dx,dy);
        disp(dx*dy*sqrt(sum(sum(b.^2))))
        err(i) = dx*dy*sqrt(sum(sum(b.^2)));
        figure(1);
        surf(xx,yy,u,'edgecolor','none');
        error = abs(u-u_exact(xx,yy));
        figure(2); clf;
        surf(xx,yy,error,'edgecolor','none');
        pause(0.1);
    end
    figure(1);
    surf(xx,yy,u,'edgecolor','none');
    error = abs(u-u_exact(xx,yy));
    figure(2); clf;
    surf(xx,yy,error,'edgecolor','none');
    L1_norm(i) = dx*dy*sum(sum(error));
    L2_norm(i) = dx*dy*sum(sum(error.^2));
    max_norm(i) = max(max(error));
    i = i+1;
    %pause()
end
figure(3); clf;
loglog(1./NN,L1_norm)
