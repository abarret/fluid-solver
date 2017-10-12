
function vals = projectSolve(vals)
global Nx; global Ny;
global rho; global dt;
global mu; global Lx; global Ly;
global dx; global dy;
fu = vals(1:(Nx)*Ny);
fv = vals(Nx*Ny+(1:Nx*Ny));
fp = vals(Nx*Ny+Nx*Ny+(1:Nx*Ny));
fu = reshape(fu, Ny, Nx);
fv = reshape(fv, Ny, Nx);
fp = reshape(fp, Ny, Nx);
% First solve for velocity
u = gmres(@velSolveU,fu(:),[],1.0e-2,25);
u = reshape(u,Ny,Nx);
u = [u, u(:,1)];
v = gmres(@velSolveV,fv(:),[],1.0e-2,25);
v = reshape(v,Ny,Nx);
v = [v; v(1,:)];
% Project velocity onto certain field
rhs_p = -rho/dt*(fp+DivergenceStoC(u,v,dx,dy));
% Solve poisson problem for phi
L = ichol(sparse(getA(@poissonSolve)));
%phi = gmres(@poissonSolve,rhs_p(:),1,1.0e-8,1000);
phi = pcg(@poissonSolve,rhs_p(:),1.0e-2,25,L,L');
%[phi, tol, iter] = PoissonSolveSOR(zeros(Nx,Ny),rhs_p,500,1.0e-8,1.25);
%fprintf('\n\n PoissonSolveSOR finished after %d iterations with final rel tol %f \n\n', iter, tol);
%phi = gmres(@poissonSolve,rhs_p(:),50,1.0e-2,500,L,L');
phi = reshape(phi,Ny,Nx);
% Take appropriate gradients...
[grad_phi_x,grad_phi_y] = GradCtoS(phi,dx,dy);
u = u - dt/rho*grad_phi_x;
v = v - dt/rho*grad_phi_y;
u = u(:,1:end-1);
v = v(1:end-1,:);
p = phi - dt/rho*mu/2*LaplacianCenter(phi,dx,dy);
vals = [u(:); v(:); p(:)];
end

function A = getA(f)
global Nx; global Ny;
A = zeros(Nx*Ny,Nx*Ny);
for i = 1:Nx*Ny
    x = zeros(Nx*Ny,1); x(i) = 1.0;
    A(:,i) = f(x);
end
end

function vals = velSolveU(vals)
global Nx; global Ny;
global dt; global rho;
global mu; global dx;
global dy;
u = reshape(vals, Ny, Nx);
u = [u, u(:,1)];
vals = rho/dt*u-mu/2*LaplacianSideX(u,dx,dy);
vals = vals(:,1:end-1);
vals = vals(:);
end

function vals = velSolveV(vals)
global Nx; global Ny;
global dt; global rho;
global mu; global dx;
global dy;
v = reshape(vals, Ny, Nx);
v = [v; v(1,:)];
vals = rho/dt*v-mu/2*LaplacianSideY(v,dx,dy);
vals = vals(1:end-1,:);
vals = vals(:);
end

function vals = poissonSolve(vals)
global Nx; global Ny;
global dx; global dy;
phi = reshape(vals, Ny, Nx);
vals = -LaplacianCenter(phi,dx,dy);
vals = vals(:);
end
