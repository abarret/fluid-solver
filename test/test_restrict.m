clear; close all;
addpath('..');

f = @(x,y) sin(pi*x).*sin(pi*y);
L = 2;
N = 16; dx = L/N;
x_side = 0:dx:L;
y_side = 0:dx:L;
x = 0.5*(x_side(2:end)+x_side(1:end-1));
y = 0.5*(y_side(2:end)+y_side(1:end-1));
[xx,yy] = meshgrid(x,y);
x_coarse_s = 0:dx*2:L; y_coarse_s = 0:dx*2:L;
x_coarse = 0.5*(x_coarse_s(2:end)+x_coarse_s(1:end-1));
y_coarse = 0.5*(y_coarse_s(2:end)+y_coarse_s(1:end-1));
[xx_coarse,yy_coarse] = meshgrid(x_coarse,y_coarse);
f_vals = f(xx,yy);

f_vals_coarse = restrict(f_vals,2,0,0);
figure(1); clf; hold on;
surf(xx,yy,f_vals,'edgecolor','none');
plot3(xx_coarse,yy_coarse,f_vals_coarse,'rx');
figure(2); clf;
surf(xx_coarse,yy_coarse,f_vals_coarse,'edgecolor','none');