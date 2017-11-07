function [u,r] = vCycle(u,b,dx,dy,tol,nMax)
[r,c] = size(u);
if r == 4 || c == 4
    % Solve problem exactly
    u = Gauss_Seidel_Poisson(u,b,tol,nMax,dx,dy);
    r = b - LaplacianCenter(u,dx,dy);
else
    % Solve current problem
    u_new = Gauss_Seidel_Poisson(u,b,tol,nMax,dx,dy);
    b_new = b - LaplacianCenter(u_new,dx,dy);
    % Restrict solution
    b_new = restrict(b_new,0,dx,dy);
    [u_temp,~] = vCycle(zeros(size(b_new)),b_new,2*dx,2*dy,tol,nMax);
    u = u_new + interpolate(u_temp,0,dx,dy);
    r = b - LaplacianCenter(u,dx,dy);
    u = Gauss_Seidel_Poisson(u,b,tol,nMax,dx,dy);
end
% function u_new = vCycle(u,b,depth,dx,dy)
%     v = cell(1,3);
%     f = cell(2,3);
%     v{1} = Gauss_Seidel_Poisson(u,b,1.0e-2,10,dx,dy);
%     f{1,1} = b - LaplacianCenter(v{1},dx,dy);
%     f{2,2} = restrict(f{1,1},0,dx,dy);
%     v{2} = Gauss_Seidel_Poisson(zeros(size(f{2,2})),f{2,2},1.0e-2,10,dx*2,dy*2);
%     f{1,2} = f{2,2} - LaplacianCenter(v{2},2*dx,2*dy);
%     f{2,3} = restrict(f{1,2},0,dx,dy);
%     v{3} = Gauss_Seidel_Poisson(zeros(size(f{2,3})),f{2,3},1.0e-2,10,dx*4,dy*4);
%     % Go back up
%     v{2} = v{2} + interpolate(v{3},0,dx,dy);
%     v{2} = Gauss_Seidel_Poisson(v{2},f{2,2},1.0e-2,10,dx*2,dy*2);
%     v{1} = v{1} + interpolate(v{2},0,dx,dy);
%     u_new = v{1};
%     u_new = Gauss_Seidel_Poisson(u_new,b,1.0e-2,10,dx,dy);
% end