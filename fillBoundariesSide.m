function [u_new, v_new] = fillBoundariesSide(u, v, gcw)
% Fills ghost cells using periodic boundary conditions
% for side centered vector field (u,v)

[r,c] = size(u);
u_new = zeros(r+2*gcw,c+2*gcw);
v_new = zeros(c+2*gcw,r+2*gcw);
u_new(gcw+1:end-gcw,gcw+1:end-gcw) = u;
v_new(gcw+1:end-gcw,gcw+1:end-gcw) = v;
% Fill bottom rows
u_new(1:gcw,gcw+1:end-gcw) = u_new(end-2*gcw+1:end-gcw,gcw+1:end-gcw);
v_new(1:gcw,gcw+1:end-gcw) = v_new(end-2*gcw:end-gcw-1,gcw+1:end-gcw);
% Fill left rows
u_new(1:end-gcw,1:gcw) = u_new(1:end-gcw,end-2*gcw:end-gcw-1);
v_new(1:end-gcw,1:gcw) = v_new(1:end-gcw,end-2*gcw+1:end-gcw);
% Fill right rows
u_new(1:end-gcw,end-gcw+1:end) = u_new(1:end-gcw,gcw+2:2*gcw+1);
v_new(1:end-gcw,end-gcw+1:end) = v_new(1:end-gcw,gcw+1:2*gcw);
% Fill top rows
u_new(end-gcw+1:end,:) = u_new(gcw+1:2*gcw,:);
v_new(end-gcw+1:end,:) = v_new(gcw+2:2*gcw+1,:);
end
