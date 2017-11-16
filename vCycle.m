function [u,r] = vCycle(u,b,dx,dy,tol,nMax)
% Does a vCycle on u for Poisson problem. Uses
% Gauss Seidel as smoother.
% INPUTS:
%   u    : Initial guess
%   b    : RHS vector
%   dx   : grid spacing in x-direction
%   dy   : grid spacing in y-direction
%   tol  : tolerance for smoother
%   nMax : max iterations for smoother
% OUTPUTS:
%   u    : final solution
%   r    : residual
[r,c] = size(u);
if r == 4 || c == 4
    % Solve problem exactly
    u = Gauss_Seidel_Poisson(u,b,1.0e-12,500,dx,dy);
    r = b - LaplacianCenter(u,dx,dy);
else
    % Smooth on current level.
    u_new = Gauss_Seidel_Poisson(u,b,0.0,nMax,dx,dy);
    b_new = b - LaplacianCenter(u_new,dx,dy);
    % Restrict solution
    b_new = restrict(b_new,0,dx,dy);
    [u_temp,~] = vCycle(zeros(size(b_new)),b_new,2*dx,2*dy,tol,nMax);
    u = u_new + interpolate(u_temp,0,dx,dy);
    r = b - LaplacianCenter(u,dx,dy);
    u = Gauss_Seidel_Poisson(u,b,0.0,nMax,dx,dy);
end
end
