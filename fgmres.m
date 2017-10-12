function [x, r0_norm, out, in] = fgmres(A,b,x0,M,rtol,nMax_in,nMax_out)
% Solves A*x = b using the preconditioner M using a flexible gmres
% algorithm.
% Inputs:
%   A        : function handle that computes A*x
%   b        : rhs vector
%   x0       : initial guess for solution.
%   M        : preconditioner that computes M^-1*b (NOTE preconditioner 
%              acts on the residual, not the solution)
%   rtol     : tolerance for solver
%   nMax_in  : Max number of inner solver iterations.
%   nMax_out : Max number of outer solver interations.
for out = 1:nMax_out
    % initialize vectors
    r0 = b-A(x0);
    beta = norm(r0,2);
    if(beta < rtol)
        % initial guess is within tolerance
        r0_norm = beta;
        x = x0;
        in = 0;
        return;
    end
    v = zeros(length(x0),nMax_in+1);
    v(:,1) = r0/beta;
    Z = zeros(length(x0),nMax_in);
    H = zeros(nMax_in+1,nMax_in);
    for in = 1:nMax_in
        Z(:,in) = M(v(:,in));
        w = A(Z(:,in));
        for i = 1:in
            H(i,in) = w'*v(:,i);
            w = w - H(i,in)*v(:,i);
        end
        H(in+1,in) = norm(w,2);
        v(:,in+1) = w/H(in+1,in);
        % Solve least squares...
        e1 = zeros(in+1,1); e1(1) = 1;
        y = H(1:in+1,1:in)\(beta*e1);
        x = x0 + Z(:,1:in)*y;
        % Compute residual
        r0 = b-A(x);
        r0_norm = norm(r0,2);
%        fprintf('\n\n fgmres at inner iteration %d and outer iteration %d \n current residual norm %f\n\n', in, out, r0_norm);
        if (r0_norm < rtol)
            return;
        end
    end
    % Hasn't converged in inner iteration. Update states for next outer
    % iteration
    x0 = x;
end