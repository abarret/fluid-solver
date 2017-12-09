clear;
% We need a few global variables. These should be double checked
% for consistency. Some are inputs to functions...
global Lx; global Ly;
global dx; global dy;
global rho;
global mu; global mup;
global lambda;
global Nx; global Ny;
global dt; global is_periodic;
counter = 0;
% Flag to check periodicity. Not currently used.
is_periodic = 1;
% Convergence tests on N....
for N = [32 64 128]
counter = counter + 1;
Nx = N; Ny = N;
Lx = 1; Ly = 1;
dt = Lx/(2*Nx); T = 0.5;
rho = 1.0;%1.0e-8;
mu = 1.0e-2; mup = 0.0;
lambda = 0.1;

% Taylor vortices as exact solution.
u_exact = @(x,y,t) 1-2*exp(-8*pi*pi*mu/rho*t)*sin(2*pi*(y-t)).*cos(2*pi*(x-t));
v_exact = @(x,y,t) 1+2*exp(-8*pi*pi*mu/rho*t)*cos(2*pi*(y-t)).*sin(2*pi*(x-t));
p_exact = @(x,y,t) -exp(-16*pi*pi*mu/rho*t)*(cos(4*pi*(x-t))+cos(4*pi*(y-t)));

% Force functions should be split.
%f = @(x,y,t) {2*sin(x).*cos(y); -2*cos(x).*sin(y)};
fx = @(x,y,t) 0*x;
fy = @(x,y,t) 0*x;

dx = Lx/Nx; dy = Ly/Ny;
x_side = 0:dx:Lx; y_side = 0:dy:Ly;
x_cent = 0.5*(x_side(2:end)+x_side(1:end-1));
y_cent = 0.5*(y_side(2:end)+y_side(1:end-1));
[xx_cent,yy_cent] = meshgrid(x_cent,y_cent);
[xx_sidex,yy_sidex] = meshgrid(x_side,y_cent);
[xx_sidey,yy_sidey] = meshgrid(x_cent,y_side);
ilower = [0 0];
iupper = [Nx Ny];
xlow = [0 0];
xup = [Lx Ly];
u = u_exact(xx_sidex,yy_sidex,0); v = v_exact(xx_sidey,yy_sidey,0);
p = p_exact(xx_cent,yy_cent,0);
t = 0; iter = 1;
draw_freq = 10;
draw_stuff(u, v, p, t);
pause(0.01);
while(abs(T-t) > 1.0e-10)
    % Calculate next time step. This should be done using a CFL condition.
    % Write a function called getTimeStep()?
    dt = min(dt, T-t);
    % Advance solution in time
    [u, v, p] = advance_in_time(u,v,p,fx,fy,t,dt);
    t = t + dt
    if((mod(iter, draw_freq) == 0) || (abs(T-t) <= 1.0e-10))
        draw_stuff(u, v, p, t);
        pause(0.01);
    end
    iter = iter+1;
end
u_exact_vals = u_exact(xx_sidex,yy_sidex,t);
v_exact_vals = v_exact(xx_sidey,yy_sidey,t);
p_exact_vals = p_exact(xx_cent,yy_cent,t-dt/2);
error_u = abs(u_exact_vals-u);
error_v = abs(v_exact_vals-v);
error_p = abs(p_exact_vals-p);
L1_u(counter) = dx*dy*sum(sum(error_u));
L2_u(counter) = dx*dy*sqrt(sum(sum(error_u.^2)));
max_u(counter) = max(max(error_u));
L1_v(counter) = dx*dy*sum(sum(error_v));
L2_v(counter) = dx*dy*sqrt(sum(sum(error_v.^2)));
max_v(counter) = max(max(error_v));
L1_p(counter) = dx*dy*sum(sum(error_p));
L2_p(counter) = dx*dy*sqrt(sum(sum(error_p.^2)));
max_p(counter) = max(max(error_p));
end

% External functions
% I think this is actually in a seperate function.
% Probably can be deleted...
function draw_stuff(u, v, p, t)
% Currently blank...
global dx; global dy;
global Lx; global Ly;
x_side = 0:dx:Lx; y_side = 0:dy:Ly;
x_cent = 0.5*(x_side(2:end)+x_side(1:end-1));
y_cent = 0.5*(y_side(2:end)+y_side(1:end-1));
[xx,yy] = meshgrid(x_cent,y_cent);
uu = sideToCell(u,v);
u_c = uu(:,:,1); v_c = uu(:,:,2);
figure(1); clf;
subplot(2,2,1);
quiver(xx,yy,u_c,v_c);
xlabel('x'); ylabel('y');
title(['velocity vector plot t = ' num2str(t)]);
subplot(2,2,2);
h = pcolor(xx,yy,p); colorbar;
set(h,'edgecolor','none');
xlabel('x'); ylabel('y');
title(['pressure plot t = ' num2str(t)]);
subplot(2,2,3);
h=pcolor(xx,yy,sqrt(u_c.^2+v_c.^2)); colorbar
set(h,'edgecolor','none');
xlabel('x'); ylabel('y');
title(['velocity mag plot t = ' num2str(t)]);
subplot(2,2,4);
h=pcolor(xx,yy,(DivergenceStoC(u,v,dx,dy)));
set(h,'edgecolor','none');
xlabel('x'); ylabel('y');
title(['div(u) plot t = ' num2str(t)]);
colorbar;
end
