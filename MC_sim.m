%% Function for running Monte Carlo simulations, given an estimated constraint set 
% Take the set. Run. Estimate failure probability and average cost of successful iterations 
% Monimoy Bujarbaruah

function [prob_fail, cost_average] = MC_sim(Xhat, iC)

    %% Define system parameters 
    [A,B,Cold,Dold,bold,Xold,Cnew,Dnew,bnew,Xnew,U,W, x_0, Q,R,N, x_ref,simsteps, nx, nu] = sys_load();
                                                       
    %% Arrays for iteration data 
    x_cl = zeros(nx,simsteps+1,iC);
    u_cl = zeros(nu,simsteps,iC);
    cost_iter = zeros(1,iC); 

    %% Things needed with estimated constraint set
    Chat_apnd = Xhat.A; 
    bhat_apnd = Xhat.b; 
    Chat = [Cold(1:4,:); Chat_apnd; Cold(5:6,:)]; 
    Dhat = [[0; 0; 0; 0]; zeros(size(Chat_apnd,1),1); [1; -1]]; 
    bhat = [bold(1:4); bhat_apnd; bold(5:6)]; 
    
    % forming the terminal set. Non-empty for sure, if used after
    % successful convergence of the main.m file
    [Xn,Pinf] = term_setRob(A,B,Chat,Dhat,bhat,Q,R,U,W,nx,nu);

    % forming matrices needed for MPC
    [matF, matG, matH, mat_c, dim_t] = obtain_matR(A,B,Chat,Dhat,bhat,Xn,N);

    %% Big iteration loop starts here
    succ_iter = 1;

    for iter_count = 1:iC
        %% Getting Started with Closed Loop MPC Problem 
        x_cl(:,1,iter_count) = x_0;  
        options = sdpsettings('solver','quadprog','verbose',0);

        %% Run an RMPC Along the Iteration (This will be feasible now) 
        % solve MPC problem, updating closed-loop information
        for i = 1:simsteps
            % solve optimization problem to minimize cost
            [M, v, ~] = solveMatrices(A, B, N, W, dim_t, matF, mat_c, matH, matG, x_cl(:,i,iter_count), x_ref, Q, R, Pinf, 1);
            % apply closed-loop input
            u_t = value(v(1));
            u_cl(:,i,iter_count) = u_t;    
            % simulate a realized disturbance 
            w_t = 2*W.b(1)*rand(2,1)-W.b(1);
            % update dynamics
            x_cl(:,i+1,iter_count) = A*x_cl(:,i,iter_count)+B*u_t + w_t;
        end

        % calculate the cost associated with the closed-loop along the iteration 
        cost_iter(1,iter_count) = calcCost(x_cl(:,:,iter_count), u_cl(:,:,iter_count), Q, R, Pinf, x_ref);


        %% Check successful iterations and then only do an average cost  
        lab_iter = simul_con(x_cl(:,:,iter_count), Xnew);
        if prod(lab_iter) == 1                                                            % successful iteration 
            cost_succIter(succ_iter) = cost_iter(1,iter_count);
            succ_iter = succ_iter + 1; 
        end

        iter_count

    end
    
    prob_fail = 1 - ((succ_iter-1)/iC);
    cost_average = mean(cost_succIter);

end

%% Results
% The test should satisfy all the expected guarantees of the estimated
% constraint set, obtained after the successful convergence of the main file

