function draw_stuff(u, v, p, sxx, syy, sxy, t)
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
% subplot(2,3,5);
% pcolor(xx,yy,sxx+syy); colorbar
% xlabel('x'); ylabel('y');
% title(['tr(S) plot t = ' num2str(t)]);
% subplot(2,3,6);
% pcolor(xx,yy,sxy); colorbar
% xlabel('x'); ylabel('y');
% title(['sxy plot t = ' num2str(t)]);
end
