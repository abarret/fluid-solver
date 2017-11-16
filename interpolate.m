function u_new = interpolate(u,ratio,dx,dy)
% Interpolates u onto finer grid using linear interpolation
% in each dimension. Used for multigrid.
% NOTE: ratio currently has no effect. Assumes we are
% interpolating onto a grid twice as large.

% Assume ratio is 2 for now
% Assume u is cell centered
[r,c] = size(u);
u_new = zeros(r*2,c*2);
u_temp = zeros(r*2,c+2);
u = fillBoundariesCenter(u,1);
% Interpolate in one dimension...
for i = 1:length(u_temp(1,:))
    u_temp(:,i) = interpolate1D(u(:,i),ratio,dx,dy);
end
% Interpolate in other dimension...
for i = 1:length(u_new(:,1))
    u_new(i,:) = interpolate1D(u_temp(i,:),ratio,dx,dy);
end
end

function u_new = interpolate1D(u,ratio,dx,dy)
% Assume ratio is 2 for now... Better interpolations can come later
% Assume u is cell centered for now....
% This only works for one spatial dimension
u_new = zeros(1,length(u(2:end-1))*2);
%u = [u(end) u u(1)];
u_new(2:2:end) = 0.25*u(3:end)+0.75*u(2:end-1);
u_new(1:2:end) = 0.75*u(2:end-1)+0.25*u(1:end-2);
end
