%% Solve the MPC Problem
% Monimoy Bujarbaruah and Charlott Vallon

function [M, v, feas] = solveMatrices(A, B, N, W, dim_t, matF, mat_c, matH, matG, x_0, x_ref, Q, R, Pinf, costflag)


    % if costflag == 1, we minimize cost. elseif costflag == 0, we only verify feasibility.

    Hs=[]; hs =[];
    for k = 1:N
        polS = W;
        Hs_ol = polS.A;  hs_ol = polS.b;
        Hs = blkdiag(Hs, Hs_ol); hs = [hs; hs_ol];
    end
    dim_a = size(Hs,1); 
    nu = size(R,2);
    nx = size(Q,2);

    % define optimization variables
    M = sdpvar(nu*N,nx*N);
    v = sdpvar(nu*N,1);
    Z = sdpvar(dim_a,dim_t);
    x = sdpvar(nx,N+1);

    for j=1:nu*N
        for k = 2*j-1:nx*N
            M(j,k) =0;
        end
    end

    % constraints for M and v
    constraints = [matF*v + Z'*hs <= mat_c + matH*x_0];
    constraints = [constraints; Z>=0];
    constraints = [constraints, matF*M + matG == Z'*Hs];
    
    % constraints for dynamics and cost
    cost = 0;
    constraints = [constraints, x(:,1) == x_0];
    for j = 1 : N
        constraints = [constraints, x(:,j+1) == A*x(:,j) + B*v(j)];
        cost = cost + (x(:,j)-x_ref)'*Q*(x(:,j)-x_ref) + v(j)*R*v(j);
    end
    cost = cost + (x(:,N+1)-x_ref)'*Pinf*(x(:,N+1)-x_ref);              

    % specify solver
    options = sdpsettings('solver','gurobi','verbose',0);

    % solve problem and check for feasibility
    if costflag == 1
        diagn=solvesdp(constraints, cost, options); %diagn.problem
    else
        diagn=solvesdp(constraints, [], options);
    end
    
    feas = diagn.problem;

end