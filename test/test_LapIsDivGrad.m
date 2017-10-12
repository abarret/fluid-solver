addpath('..');
u = @(x,y) exp(-100*((x-0.5).^2+(y-0.5).^2));%cos(2*pi*x).*sin(2*pi*y);

Lx = 1.0; Ly = 1.0;
NN = [32 64 128 256 512 1024]

L1 = 0*NN; L2 = 0*NN; max_e = 0*NN;
i = 1;
for N = NN
x_side = linspace(0,1,N+1); y_side = linspace(0,1,N+1);
x_cent = 0.5*(x_side(2:end)+x_side(1:end-1));
y_cent = 0.5*(y_side(2:end)+y_side(1:end-1));
[xx,yy] = meshgrid(x_cent,y_cent);
lap_u_approx = Laplacian(u(xx,yy), x_cent(2)-x_cent(1), y_cent(2)-y_cent(1));
[gradx,grady] = GradCtoS(u(xx,yy), x_cent(2)-x_cent(1), y_cent(2)-y_cent(1));
lap_u_gd = DivergenceStoC(gradx,grady,x_cent(2)-x_cent(1),y_cent(2)-y_cent(1));

figure(1); clf;
surf(xx,yy,lap_u_approx,'edgecolor','none'); xlabel('x'); ylabel('y');
error = abs(lap_u_approx-lap_u_gd);
figure(2); clf;
surf(xx,yy,error,'edgecolor','none'); xlabel('x'); ylabel('y');
dx = x_cent(2) - x_cent(1); dy = y_cent(2)-y_cent(1);
L1(i) = dx*dy*sum(sum(error));
L2(i) = dx*dy*sum(sum(error.^2));
max_e(i) = max(max(error));
i = i+1;
end