addpath('..');
clear;
q_fcn = @(x,y,t) heaviside(x-t-0.55).*heaviside(y-t-0.55)-heaviside(x-t-0.85).*heaviside(y-t-0.55)-heaviside(x-t-0.55).*heaviside(y-t-0.85)+heaviside(x-t-0.85).*heaviside(y-t-0.85) + exp(-100*((x-t-0.2).^2+(y-t-0.2).^2));%exp(-100*((x-t-0.5).^2+(y-t-0.5).^2));
u_fcn = @(x,y) ones(size(x));%y-0.5;
v_fcn = @(x,y) ones(size(x));%-(x-0.5);

L_x = 1.0; L_y = 1.0;
global dx; global dy;
count = 1;
for N = [32 64 128 256]
    dx = L_x/N; dy = L_y/N;
    x_side = 0:dx:L_x; y_side = 0:dy:L_y;
    x_cent = 0.5*(x_side(1:end-1)+x_side(2:end));
    y_cent = 0.5*(y_side(1:end-1)+y_side(2:end));

    [xx_cent,yy_cent] = meshgrid(x_cent,y_cent);
    [xx_sidex,yy_sidex] = meshgrid(x_side,y_cent);
    [xx_sidey,yy_sidey] = meshgrid(x_cent,y_side);

    u = u_fcn(xx_sidex,yy_sidex);
    v = v_fcn(xx_sidey,yy_sidey);
    q = q_fcn(xx_cent,yy_cent,0);

    T = 1; dt = 0.25*dx;
    t = 0; iter = 1;
    draw_freq = 10;
    draw_stuff(xx_cent,yy_cent,q,u_fcn,v_fcn,t);
    while(abs(T-t) > 1.0e-12)
        dt = min(dt, T-t);
        rhs = advect(q,u,v);
        q_temp = q+0.5*dt*rhs;
        rhs = advect(q_temp,u,v);
        q = q + dt*rhs;
        t = t + dt
        iter = iter+1;
        if((mod(iter, draw_freq) == 0) || (abs(T-t) <= 1.0e-10))
            draw_stuff(xx_cent,yy_cent,q,u_fcn,v_fcn,t);
        end
    end
    q_exact = q_fcn(xx_cent,yy_cent,t);
    err = abs(q_exact - q);
    L1(count) = dx*dy*sum(sum(err));
    L2(count) = dx*dy*sqrt(sum(sum(err.^2)));
    max_e(count) = max(max(err));
    count = count + 1;
end

function draw_stuff(x,y,q,u_fcn,v_fcn,t);
    figure(1); clf;
    subplot(1,2,1);
    h = surf(x,y,q);
    set(h,'edgecolor','none');
    subplot(1,2,2);
    quiver(x,y,u_fcn(x,y),v_fcn(x,y));
    pause(0.05);
end