function lap = LaplacianSideX(w,dx,dy)
% Computes the laplcian of side centered w (on the x axis)
gcw = 2;
% w = fillBoundaries(w,gcw);
[w, ~] = fillBoundariesSide(w,w',gcw);
lap = (w(gcw+1:end-gcw,gcw+2:end-gcw+1)-2*w(gcw+1:end-gcw,gcw+1:end-gcw)+w(gcw+1:end-gcw,gcw:end-gcw-1))/dx^2 + ...
      (w(gcw+2:end-gcw+1,gcw+1:end-gcw)-2*w(gcw+1:end-gcw,gcw+1:end-gcw)+w(gcw:end-gcw-1,gcw+1:end-gcw))/dy^2;
end