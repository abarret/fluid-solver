function [u,r] = fmg(u,b,tol,nMax,dx,dy)
% Solves the Poisson problem using a full multigrid cycle
% Uses Gauss-Seidel as the smoother
% INPUTS:
%   u     : Initial guess for solver
%   b     : RHS vector
%   tol   : tolerance level
%   nMax  : Max iterations for smoother
%   dx    : grid spacing in x-direction
%   dy    : grid spacing in y-direction
% OUTPUTS:
%   u     : Final solution
%   r     : Residual
u_temp{1} = u; b_temp{1} = b;
dx_temp = dx; dy_temp = dy;
[r,c] = size(u_temp);
i = 1;
while(r ~= 4 || c ~= 4)
    u_temp{i+1} = restrict(u_temp{i},0,dx,dy);
    b_temp{i+1} = restrict(b_temp{i},0,dx,dy);
    dx_temp = dx_temp*2;
    dy_temp = dy_temp*2;
    [r,c] = size(u_temp{i+1});
    i = i+1;
end
u_temp{end} = Gauss_Seidel_Poisson(u_temp{end},b_temp{end},tol,nMax,dx_temp,dy_temp);
for j = i-1:-1:1
    dx_temp = dx_temp/2; dy_temp = dy_temp/2;
    [u_temp{j}] = vCycle(interpolate(u_temp{j+1},0,dx,dy),b_temp{j},dx_temp,dy_temp,tol,nMax);
end
u = u_temp{1};
r = b - LaplacianCenter(u,dx,dy);
end
