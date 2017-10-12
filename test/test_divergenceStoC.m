addpath('..');
u = @(x,y) sin(2*pi*x).*sin(2*pi*y);
v = @(x,y) cos(2*pi*y).*sin(2*pi*x);
NN = [16 32 64 128 256 512 1024 2048 4096];
L1 = 0*NN; L2 = 0*NN; max_e = 0*NN;
i = 1;
for N = NN
x_side = linspace(0,1,N+1); y_side = linspace(0,1,N+1);
x_cent = 0.5*(x_side(2:end)+x_side(1:end-1));
y_cent = 0.5*(y_side(2:end)+y_side(1:end-1));
[xx_xs,yy_xs] = meshgrid(x_side,y_cent);
[xx_ys,yy_ys] = meshgrid(x_cent,y_side);
[xx, yy] = meshgrid(x_cent,y_cent);
div = DivergenceStoC(u(xx_xs,yy_xs), v(xx_ys,yy_ys), x_side(2)-x_side(1), y_side(2)-y_side(1));

%div_exact = @(x,y) 2*pi*cos(2*pi*x)-2*pi*sin(2*pi*y);
div_exact = @(x,y) 2*pi*(cos(2*pi*x).*sin(2*pi*y)-sin(2*pi*x).*sin(2*pi*y));
div_exact_vals = div_exact(xx,yy);
figure(2); clf;
surf(xx,yy,div_exact_vals-div);
error = abs(div-div_exact_vals);
dx = x_side(2) - x_side(1); dy = y_side(2)-y_side(1);
L1(i) = dx*dy*sum(sum(error));
L2(i) = dx*dy*sum(sum(error.^2));
max_e(i) = max(max(error));
i = i+1;
end
figure(1); clf; hold on;
plot(log(1./NN),log(L1),'b','linewidth',2);
a1=polyfit(log(1./NN),log(L1),1);
plot(log(1./NN),log(L2),'g','linewidth',2);
a2=polyfit(log(1./NN),log(L2),1);
plot(log(1./NN),log(max_e),'r','linewidth',2);
am=polyfit(log(1./NN),log(max_e),1);
disp(['slope for L1 = ' num2str(a1(1))]);
disp(['slope for L2 = ' num2str(a2(1))]);
disp(['slope for max = ' num2str(am(1))]);
hold off;