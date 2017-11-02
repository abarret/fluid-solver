function u_new = restrict(u,ratio,dx,dy)
[r,c] = size(u);
u_new = zeros(r/2,c/2);
u_temp = zeros(r/2,c+2);
u = fillBoundariesCenter(u,1);
% Interpolate in one dimension...
for i = 1:length(u_temp(1,:))
    u_temp(:,i) = restrict1D(u(:,i),ratio,dx);
end
% Interpolate in other dimension...
for i = 1:length(u_new(:,1))
    u_new(i,:) = restrict1D(u_temp(i,:),ratio,dx);
end
end

function u_new = restrict1D(u,ratio,dx)
u_new = 1/8*u(1:2:end-3)+3/8*u(2:2:end-2)+3/8*u(3:2:end-1)+1/8*u(4:2:end);
end