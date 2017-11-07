clear; close all;
addpath('..');

f = @(x,y) sin(pi*x).*sin(pi*y);
L = 2; i = 1;
NN = [32 64 128 256 512 1024];
for N = NN
    dx = L/N;
    x_side = 0:dx:L;
    y_side = 0:dx:L;
    x = 0.5*(x_side(2:end)+x_side(1:end-1));
    y = 0.5*(y_side(2:end)+y_side(1:end-1));
    [xx,yy] = meshgrid(x,y);
    x_fine_s = 0:dx/2:L; y_fine_s = 0:dx/2:L;
    x_fine = 0.5*(x_fine_s(2:end)+x_fine_s(1:end-1));
    y_fine = 0.5*(y_fine_s(2:end)+y_fine_s(1:end-1));
    [xx_fine,yy_fine] = meshgrid(x_fine,y_fine);
    f_vals = f(xx,yy);

    f_vals_fine = interpolate(f_vals,2,0,0);
    figure(1); clf;
    surf(xx,yy,f_vals,'edgecolor','none');
    figure(2); clf;
    surf(xx_fine,yy_fine,f_vals_fine,'edgecolor','none');
    figure(3); clf;
    surf(xx_fine,yy_fine,abs(f_vals_fine-f(xx_fine,yy_fine)),'edgecolor','none');
    L1(i) = dx*dx*sum(sum(abs(f_vals_fine-f(xx_fine,yy_fine))));
    L2(i) = dx*dx*sqrt(sum(sum((f_vals_fine-f(xx_fine,yy_fine)).^2)));
    i = i+1
end
figure(4); clf;
loglog(NN,L1);
hold on;
loglog(NN,L2,'r');