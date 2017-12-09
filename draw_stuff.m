function draw_stuff(u, v, p, t)
% Draws a 4 subplots on a single figure. Currently draws the velocity
% vector field, the magnitude of the velocity, the pressure, and the
% divergence of the velocity field.
% Requires an extra stress field that is currently not used. This is a
% relic of an oldroyd-B solve code that might come back.
% Inputs:
%   u   : side-centered x velocity
%   v   : side-centered y velocity
%   p   : cell-centered pressure
%   t   : current time
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
