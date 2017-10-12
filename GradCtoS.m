function [u, v] = GradCtoS(p, dx, dy)
% Computes the gradient of a cell centered quantity onto a side centered
% vector field.
% Inputs:
%   p  : Cell-centerd scalar field
%   dx : grid spacing in x-direction
%   dy : grid spacing in y-direction
p = fillBoundariesCenter(p, 1);
v = (p(2:end,2:end-1)-p(1:end-1,2:end-1))/dy;
u = (p(2:end-1,2:end)-p(2:end-1,1:end-1))/dx;
end