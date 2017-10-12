function [u,v,p] = advanceStokes(u,v,p,u_0,v_0,fx,fy,t,dt)
% Solves the navier stokes equations
% Inputs:
%   u   : side-centered x velocity component
%   v   : side-centered y velocity component
%   p   : cell-centered pressure (used as initial guess only)
%   u_0 : x velocity at beginning of timestep
%   v_0 : y velocity at beginning of timestep
%   fx  : function that computes the x component of the force
%   fy  : function that computes the y component of the force
%   t   : current time
%   dt  : time step
global rho; global mu;
global Lx; global Ly;
global Nx; global Ny;
global dx; global dy;
global is_periodic;
gcw = 1;
% Calculate nonlinear advection term
[u_advect, v_advect] = fillBoundariesSide(0.5*(u+u_0),0.5*(v+v_0),gcw);
uu_advect = sideToCell(u_advect,v_advect);
Nu = advect(0.5*(u+u_0),uu_advect(2:end-1,:,1),0.5*(v_advect(2:end-1,2:end)+v_advect(2:end-1,1:end-1)));
Nv = advect(0.5*(v+v_0),0.5*(u_advect(1:end-1,2:end-1)+u_advect(2:end,2:end-1)),uu_advect(:,2:end-1,2));
x_side = 0:dx:Lx; y_side = 0:dy:Ly;
x_cent = 0.5*(x_side(2:end)+x_side(1:end-1));
y_cent = 0.5*(y_side(2:end)+y_side(1:end-1));
[xx_sidex, yy_sidex] = meshgrid(x_side,y_cent);
[xx_sidey, yy_sidey] = meshgrid(x_cent,y_side);
% Evaluate forces at half timestep
f_valsx = fx(xx_sidex,yy_sidex,t+dt/2);
f_valsy = fy(xx_sidey,yy_sidey,t+dt/2);
% Form rhs vector for fgmres cycle
r1 = rho/dt*u_0+mu/2*LaplacianSideX(u_0,dx,dy)-rho*Nu+f_valsx;
r2 = rho/dt*v_0+mu/2*LaplacianSideY(v_0,dx,dy)-rho*Nv+f_valsy;
r3 = zeros(size(p));
% Remove extra boundary components from r1 and r2. Note that this only has to
% be done for periodic boundaries.
% ignore top row of r2 (y side values)
r2 = r2(1:end-1,:);
% ignore right column of r1 (x side values)
r1 = r1(:,1:end-1);
rhs_vec = [r1(:); r2(:); r3(:)];
% Remove extra boundary components from u and v. Note that this only has to
% be done for periodic boundaries.
% ignore top row of A2 (y side values)
v = v(1:end-1,:);
% ignore right column of A1 (x side values)
u = u(:,1:end-1);
x_vec = [u(:); v(:); p(:)];
% Solve system. Uses a flexible gmres algorithm alowing the preconditioner
% to change between steps.
[x_vec, res, out, in] = fgmres(@StokesSolve,rhs_vec,x_vec,@projectSolve,1.0e-5,200,5);
fprintf('\n\nfgmres converged in %d outer iterations with %d inner iterations with residual %f\n\n',out,in,res);
% Solve by LU (NOT ADVISED)
%x_vec = LUSolveStokes(@StokesSolve,rhs_vec);
% Seperate solution components
u = x_vec(1:(Nx)*Ny);
v = x_vec((Nx)*Ny+(1:Nx*(Ny)));
p = x_vec((Nx)*Ny+Nx*(Ny)+(1:Nx*Ny));
u = reshape(u,Ny,Nx);
v = reshape(v,Ny,Nx);
p = reshape(p,Ny,Nx);
% Add in boundary components for u and v. Note this only needs to be done
% for periodic boundaries.
% Add in top row for v
v = [v; v(1,:)];
% Add in right column for u
u = [u, u(:,1)];
% Normalize pressure. Note this only has to be done for boundary conditions
% that don't supply a pressure value (i.e. NOT traction boundaries)
p = p - dx*dy*sum(p);
end

function x = LUSolve(f, b)
% Takes a function f that computes matrix vector product and computes the
% corresponding matrix explicitly. Then solves system via linsolve.
A = zeros(length(b));
for i = 1:length(b)
    x = zeros(size(b)); x(i) = 1.0;
    A(:,i) = f(x);
end
x = linsolve(A,b);
end

function vals = StokesSolve(vals)
% Computes matrix vector product for the Stokes solve iteration
% Inputs:
%   vals : rhs vector stacked with [u(:); v(:); p(:)]
% Outputs:
%   vals : A*[u(:); v(:); p(:)] with
%              A = [rho*dt*I-mu/2*Lap^x          0          Grad^x;
%                        0              rho*dt*I-mu/2*Lap^y Grad^y;
%                      -Div^x                 -Div^y           0]
global Nx; global Ny;
global rho; global dt;
global mu; global Lx; global Ly;
global dx; global dy;
% u = vals(1:(Nx+1)*Ny);
% v = vals((Nx+1)*Ny+(1:Nx*(Ny+1)));
% p = vals((Nx+1)*Ny+Nx*(Ny+1)+(1:Nx*Ny));
% u = reshape(u,Ny,Nx+1);
% v = reshape(v,Ny+1,Nx);
% p = reshape(p,Ny,Nx);
% We removed certain boundary values for periodic boundary conditions.
% These need to be added in ONLY in this case.
u = vals(1:(Nx)*Ny);
v = vals((Nx)*Ny+(1:Nx*(Ny)));
p = vals((Nx)*Ny+Nx*(Ny)+(1:Nx*Ny));
u = reshape(u,Ny,Nx);
v = reshape(v,Ny,Nx);
p = reshape(p,Ny,Nx);
% Add in top row for v
v = [v; v(1,:)];
% Add in right column for u
u = [u, u(:,1)];
[gradPx,gradPy] = GradCtoS(p,dx,dy);
A1 = rho/dt*u-mu/2*LaplacianSideX(u,dx,dy)+gradPx;
A2 = rho/dt*v-mu/2*LaplacianSideY(v,dx,dy)+gradPy;
A3 = -DivergenceStoC(u,v,dx,dy);
% Remove certain boundary components. Note this only needs to be done with
% periodic boundary conditions.
% ignore top row of A2 (y side values)
A2 = A2(1:end-1,:);
% ignore right column of A1 (x side values)
A1 = A1(:,1:end-1);
vals = [A1(:); A2(:); A3(:)];
end