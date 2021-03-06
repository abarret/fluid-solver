function lap = LaplacianCenter(w,dx,dy)
% Computes the Laplacian of a cell centered quantity using second order
% finite differences.
% Inputs:
%   w  : cell centered scalar valued function.
%   dx : grid spacing in x-direciton.
%   dy : grid spacing in y-direction.
gcw = 2;
w = fillBoundariesCenter(w,gcw);
lap = (w(gcw+1:end-gcw,gcw+2:end-gcw+1)-2*w(gcw+1:end-gcw,gcw+1:end-gcw)+w(gcw+1:end-gcw,gcw:end-gcw-1))/dx^2 + ...
      (w(gcw+2:end-gcw+1,gcw+1:end-gcw)-2*w(gcw+1:end-gcw,gcw+1:end-gcw)+w(gcw:end-gcw-1,gcw+1:end-gcw))/dy^2;
end