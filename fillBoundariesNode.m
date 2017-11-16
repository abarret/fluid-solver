function u_new = fillBoundariesNode(u,gcw)
% Fills boundaries for node centered data.
% NOTE: Assumes periodicity. Top and right rows should be removed PRIOR to
% calling this function.
[r,c] = size(u);
u_new = zeros(r+2*gcw,c+2*gcw);
u_new(gcw+1:end-gcw,gcw+1:end-gcw) = u;
% Fill top rows
u_new(1:gcw,gcw+1:end-gcw) = u_new(end-2*gcw+1:end-gcw,gcw+1:end-gcw);
% Fill left columns
u_new(1:end-gcw,1:gcw) = u_new(1:end-gcw,end-2*gcw+1:end-gcw);
% Fill right columns
u_new(1:end-gcw,end-gcw+1:end) = u_new(1:end-gcw,gcw+1:2*gcw);
% Fill bottom rows
u_new(end-gcw+1:end,:) = u_new(gcw+1:2*gcw,:);
end