function vals = projectSolve(vals)
% A projection preconditioner for the incompressible stokes equations
% solved via a projection method. Note that in a fgmres solver, the
% preconditioner acts on the residual vector.
% Inputs:
%   vals : residual vector.
% Outputs:
%   vals : M^-1*vals where M is a projection preconditioner, i.e. we solve
%          Stokes equations without pressure, then project velocity on a
%          divergence free space.
global Nx; global Ny;
global rho; global dt;
global mu; global Lx; global Ly;
global dx; global dy;
% reshape RHS for the projection solver. This is probably not necessary...
fu = vals(1:(Nx)*Ny);
fv = vals(Nx*Ny+(1:Nx*Ny));
fp = vals(Nx*Ny+Nx*Ny+(1:Nx*Ny));
fu = reshape(fu, Ny, Nx);
fv = reshape(fv, Ny, Nx);
fp = reshape(fp, Ny, Nx);
% First solve for velocity
%u = gmres(@velSolveU,fu(:),[],1.0e-2,25);
[u,~,res,iter] = pcg(@velSolveU,fu(:),1.0e-2,25,[]);
u = reshape(u,Ny,Nx);
fprintf("Solved velocity equation in x-direction with residual %f in %d iterations.\n",res,iter);
% Add in extra row for boundary conditions. Note this is only needed in
% periodic boundary condtions.
u = [u, u(:,1)];
%v = gmres(@velSolveV,fv(:),[],1.0e-2,25);
[v,~,res,iter] = pcg(@velSolveV,fv(:),1.0e-2,25,[]);
v = reshape(v,Ny,Nx);
fprintf("Solved velocity equation in y-direction with residual %f in %d iterations.\n",res,iter);
% Add in extra column for boundary conditions. Note this is only needed in
% periodic boundary condtions.
v = [v; v(1,:)];
% Project velocity onto certain field. Note that we don't require the
% divergence of the velocity field to be 0, but whatever is in fp
rhs_p = -rho/dt*(fp+DivergenceStoC(u,v,dx,dy));
% Solve poisson problem for phi
% Preconditioner....
%L = ichol(sparse(getA(@poissonSolve)));
% pcg/gmres doesn't make much of a difference. We really need a better
% preconditioner here. Some kind of multigrid algorithm would be ideal.
[phi,~,res,iter] = pcg(@poissonSolve,rhs_p(:),1.0e-2,25,@fmg_wrapper);
%phi = gmres(@poissonSolve,rhs_p(:),50,1.0e-2,500,L,L');
phi = reshape(phi,Ny,Nx);
fprintf("Solved pressure equation with residual %f in %d iterations.\n",res,iter);
% Take appropriate gradients...
[grad_phi_x,grad_phi_y] = GradCtoS(phi,dx,dy);
u = u - dt/rho*grad_phi_x;
v = v - dt/rho*grad_phi_y;
u = u(:,1:end-1);
v = v(1:end-1,:);
p = phi - dt/rho*mu/2*LaplacianCenter(phi,dx,dy);
% Do we need to normalize pressure here?
vals = [u(:); v(:); p(:)];
end

function A = getA(f)
% Gets the matrix corresponding to the function f that computes A*x. This
% is used in computing the preconditioner. We really don't want to have to
% use this....
global Nx; global Ny;
A = zeros(Nx*Ny,Nx*Ny);
for i = 1:Nx*Ny
    x = zeros(Nx*Ny,1); x(i) = 1.0;
    A(:,i) = f(x);
end
end

function vals = velSolveU(vals)
% This is the velocity subdomain solver for the x direction. Perhaps we
% want to combine the two?
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
% This is the velocity subdomain solver for the y direction. Perhaps we
% want to combine the two?
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
% Solves the poisson equation for the pressure. Desperately need a good
% preconditioner for this, preferable a multigrid solver...
global Nx; global Ny;
global dx; global dy;
phi = reshape(vals, Ny, Nx);
vals = -LaplacianCenter(phi,dx,dy);
vals = vals(:);
end