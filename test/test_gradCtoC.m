addpath('..');
u = @(x,y) exp(-20*(x.^2+y.^2));
NN = [16 32 64 128 256 512 1024];
L1x = 0*NN; L2x = 0*NN; max_ex = 0*NN;
i = 1;
for N = NN
x_side = linspace(-1,1,N+1); y_side = linspace(-1,1,N+1);
x = 0.5*(x_side(2:end)+x_side(1:end-1));
y = 0.5*(y_side(2:end)+y_side(1:end-1));
[xx,yy] = meshgrid(x,y);
[gradx, grady] = Grad(u(xx,yy), x(2)-x(1), y(2)-y(1));

gradx_exact = @(x,y) -40*exp(-20*(x.^2+y.^2)).*x;
grady_exact = @(x,y) -40*exp(-20*(x.^2+y.^2)).*y;
gradx_exact_vals = gradx_exact(xx,yy);
grady_exact_vals = grady_exact(xx,yy);
figure(1); clf;
surf(xx,yy,gradx,'edgecolor','none'); xlabel('x'); ylabel('y');
errorx = abs(gradx-gradx_exact_vals);
errory = abs(grady-grady_exact_vals);
figure(2); clf;
surf(xx,yy,errorx,'edgecolor','none'); xlabel('x'); ylabel('y');
dx = x(2) - x(1); dy = y(2)-y(1);
L1x(i) = dx*dy*sum(sum(errorx));
L2x(i) = dx*dy*sum(sum(errorx.^2));
max_ex(i) = max(max(errorx));
L1y(i) = dx*dy*sum(sum(errory));
L2y(i) = dx*dy*sum(sum(errory.^2));
max_ey(i) = max(max(errory));
i = i+1;
end