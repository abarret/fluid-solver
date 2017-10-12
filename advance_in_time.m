function [u, v, p, sxx, syy, sxy] = advance_in_time(u,v,p,sxx,syy,sxy,fx,fy,t,dt)
% Solves Navier Stokes Oldroyd-B using a first order Godunov splitting scheme...
% u, v         : velocities
% sxx, syy, sxy: components of extra stress tensor
% f            : is background force
% dt           : time-step

MAX_FP_ITER=3;
% initialize solution vectors
u_k0 = u; v_k0 = v;
% Fixed pont iteration
for k = 1:MAX_FP_ITER
    % Stokes solve
    [u, v, p] = advanceStokes(u,v,p,u_k0,v_k0,fx,fy,t,dt);
end
end