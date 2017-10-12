function div = DivergenceCtoC(u,v,dx,dy)
% Computes the divergence of cell centered [u;v]
% Using second order finite differences

u = fillBoundariesCenter(u,1);
v = fillBoundariesCenter(v,1);
div = (u(2:end-1,3:end)-u(2:end-1,1:end-2))/(2.0*dx)...
    +(v(3:end,2:end-1)-v(1:end-2,2:end-1))/(2.0*dy);
end