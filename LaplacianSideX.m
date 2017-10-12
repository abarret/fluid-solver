function lap = LaplacianSideX(w,dx,dy)
% Computes the Laplacian of a side centered quantity (x-side) using second 
% order finite differences.
% Inputs:
%   w  : side centered (x-side) scalar valued function.
%   dx : grid spacing in x-direciton.
%   dy : grid spacing in y-direction.
gcw = 2;
[w, ~] = fillBoundariesSide(w,w',gcw);
lap = (w(gcw+1:end-gcw,gcw+2:end-gcw+1)-2*w(gcw+1:end-gcw,gcw+1:end-gcw)+w(gcw+1:end-gcw,gcw:end-gcw-1))/dx^2 + ...
      (w(gcw+2:end-gcw+1,gcw+1:end-gcw)-2*w(gcw+1:end-gcw,gcw+1:end-gcw)+w(gcw:end-gcw-1,gcw+1:end-gcw))/dy^2;
end