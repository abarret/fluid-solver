function div = DivergenceStoC(u,v,dx,dy)
% Computes the cell centered divergence of side centered [u;v]
% Using second order finite differences

div = (u(1:end,2:end)-u(1:end,1:end-1))/(dx)...
    +(v(2:end,1:end)-v(1:end-1,1:end))/(dy);
end