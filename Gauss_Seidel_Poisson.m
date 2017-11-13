function [u,r] = Gauss_Seidel_Poisson(u,b,tol,Nmax,dx,dy)
[r,c] = size(u);
u_temp = zeros(r,c);
u_temp_2 = zeros(r,c);
u_temp = u;
err = b - LaplacianCenter(u_temp,dx,dy);
err = dx*dy*sqrt(sum(sum(err.^2)));
n = 1;
while(err > tol)
    u_temp = fillBoundariesCenter(u_temp,1);
    u_temp_2 = fillBoundariesCenter(u_temp_2,1);
    for i = 1:r
        for j = 1:c
            % Add one to each dimension to account for ghost cells...
            u_temp_2(i+1,j+1) = (u_temp_2(i,j+1) + u_temp(i+2,j+1) + u_temp_2(i+1,j) + u_temp(i+1,j+2) - dx*dy*b(i,j)) / 4;
            %u_temp_2(i+1,j+1) = (u_temp(i,j+1) + u_temp(i+2,j+1) + u_temp(i+1,j) + u_temp(i+1,j+2) - dx*dy*b(i,j)) / 4;
        end
    end
    u_temp = u_temp_2(2:end-1,2:end-1);
    u_temp_2 = u_temp;
    err = b - LaplacianCenter(u_temp,dx,dy);
    err = dx*dy*sum(sum(abs(err)));
%    fprintf("Residual at step %d is %f \n",n,err);
    if (n >= Nmax)
        break;
    end
    n = n + 1;
end
u = u_temp;
r = b - LaplacianCenter(u,dx,dy);
end