function u_new = fillBoundariesCenter(u, gcw)
% Fills ghost cells using periodic boundary conditions

[r,c] = size(u);
u_new = zeros(r+2*gcw,c+2*gcw);
u_new(gcw+1:end-gcw,gcw+1:end-gcw) = u;
% Fill top rows
u_new(1:gcw,gcw+1:end-gcw) = u_new(end-2*gcw+1:end-gcw,gcw+1:end-gcw);
% Fill left rows
u_new(1:end-gcw,1:gcw) = u_new(1:end-gcw,end-2*gcw+1:end-gcw);
% Fill right rows
u_new(1:end-gcw,end-gcw+1:end) = u_new(1:end-gcw,gcw+1:2*gcw);
% Fill top rows
u_new(end-gcw+1:end,:) = u_new(gcw+1:2*gcw,:);
end