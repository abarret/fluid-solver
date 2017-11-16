function [u, v, p, sxx, syy, sxy] = advance_in_time(u,v,p,sxx,syy,sxy,fx,fy,t,dt)
% Solves Navier Stokes Oldroyd-B using a first order Godunov splitting scheme...
% u, v         : velocities
% p            : pressure: used as initial guess for solver
% sxx, syy, sxy: components of extra stress tensor (not currently used)
% fx, fy       : is background force functions
% t            : current time
% dt           : time-step

% Fixed point iterations. This should probably be an input somewhere.
% 2 iterations corresponds to explicit midpoint for advection and 
%              Crank-Nicolson for viscous terms.
% 1 iteration  corresponds to forward Euler for advection and
%              Crank-Nicolson for viscous terms.
% Usually use 3 iterations for improved stability
MAX_FP_ITER=3;
% initialize solution vectors
u_k0 = u; v_k0 = v;
% Fixed pont iteration
for k = 1:MAX_FP_ITER
    % Stokes solve. Note pressure is *only* used for initial guess
    % for iterative solver
    [u, v, p] = advanceStokes(u,v,p,u_k0,v_k0,fx,fy,t,dt);
end
end
