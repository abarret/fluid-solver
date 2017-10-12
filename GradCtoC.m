function [u, v] = GradCtoC(w, dx, dy)
w_new = fillBoundariesCenter(w, 1);
u = (w_new(2:end-1,3:end)-w_new(2:end-1,1:end-2))/(2.0*dx);
v = (w_new(3:end,2:end-1)-w_new(1:end-2,2:end-1))/(2.0*dy);
end