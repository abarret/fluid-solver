addpath('..');
u = @(x,y) sin(2*pi*x).*cos(2*pi*y);
NN = [16 32 64 128 256 512 1024 2048 4096];
L1x = 0*NN; L2x = 0*NN; max_ex = 0*NN;
i = 1;
for N = NN
x_side = linspace(-2,2,N+1); y_side = linspace(-1,1,N+1);
x_cent = 0.5*(x_side(2:end)+x_side(1:end-1));
y_cent = 0.5*(y_side(2:end)+y_side(1:end-1));
[xx_cent,yy_cent] = meshgrid(x_cent,y_cent);
[xx_sidex,yy_sidex] = meshgrid(x_side,y_cent);
[xx_sidey,yy_sidey] = meshgrid(x_cent,y_side);
[gradx, grady] = GradCtoS(u(xx_cent,yy_cent), x_cent(2)-x_cent(1), y_cent(2)-y_cent(1));

gradx_exact = @(x,y) 2*pi*cos(2*pi*x).*cos(2*pi*y);
grady_exact = @(x,y) -2*pi*sin(2*pi*x).*sin(2*pi*y);
gradx_exact_vals = gradx_exact(xx_sidex,yy_sidex);
grady_exact_vals = grady_exact(xx_sidey,yy_sidey);
figure(1); clf;
surf(xx_sidex,yy_sidex,gradx,'edgecolor','none'); xlabel('x'); ylabel('y');
errorx = abs(gradx-gradx_exact_vals);
errory = abs(grady-grady_exact_vals);
figure(2); clf;
surf(xx_sidex,yy_sidex,errorx,'edgecolor','none'); xlabel('x'); ylabel('y');
dx = x_cent(2) - x_cent(1); dy = y_cent(2)-y_cent(1);
L1x(i) = dx*dy*sum(sum(errorx));
L2x(i) = dx*dy*sum(sum(errorx.^2));
max_ex(i) = max(max(errorx));
L1y(i) = dx*dy*sum(sum(errory));
L2y(i) = dx*dy*sum(sum(errory.^2));
max_ey(i) = max(max(errory));
i = i+1;
end

figure(1); clf; hold on;
plot(log(1./NN),log(L1y),'b','linewidth',2);
a1=polyfit(log(1./NN),log(L1y),1);
plot(log(1./NN),log(L2y),'g','linewidth',2);
a2=polyfit(log(1./NN),log(L2y),1);
plot(log(1./NN),log(max_ey),'r','linewidth',2);
am=polyfit(log(1./NN),log(max_ey),1);
disp(['slope for L1 = ' num2str(a1(1))]);
disp(['slope for L2 = ' num2str(a2(1))]);
disp(['slope for max = ' num2str(am(1))]);
hold off;