function [u, v] = GradCtoS(p, dx, dy)
p = fillBoundariesCenter(p, 1);
v = (p(2:end,2:end-1)-p(1:end-1,2:end-1))/dy;
u = (p(2:end-1,2:end)-p(2:end-1,1:end-1))/dx;
end