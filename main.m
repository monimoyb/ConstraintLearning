%% Learning to Satisfy Unknown Constraints in Iterative MPC
% The code will stop at any one of the "keyboard" points. 
% Stopping conditions are:
% 1) if cvx hull estimate is chosen: robustness cerificate obtained or back-up set passes proposition 1 check
% 2) svm based estimate is chosen: set passes proposition 1 check
%% When to run this code?
% All the data have been saved in the .mat files in this repository, having
% run this code for the 3 cases shown in
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9303765. Run
% this code only if new data is needed. Otherwise go to file compare_code.m
%%
clc
clear all; close all;
rng(2);                                                                      % needed for comparison between different epsilons and betas 

%% Define system parameters and check for feasibility
[A,B,Cold,Dold,bold,Xold,Cnew,Dnew,bnew,Xnew,U,W, x_0, Q,R,N, x_ref,simsteps, nx, nu] = sys_load();

%% Set the epsilon and beta values here. That'll decide verification iteration count 
epsi = 0.3;
beta = 1e-7;
Nlb_probVer = ceil(log(1/beta)/log(1/(1-epsi))); 
iC = Nlb_probVer + 500;                                            % total number of iterations. Kept large, but we will terminate before

%% Some key parameters 
%%%Choose these
step_size = 0.5;                                                         % needed for set scaling to ensure MPC feasibility
cvx_flag = 0;                                                             % pick this. If 0 then SVM method. If 1 then cvx hull method. 
count_succ = 0;                                                         % successful iteration count for Proposition 1 check
%%% These are fixed
rob_flg = 0;                                                               % checks for a robustness cerificate (cvx hull method. )
prob_flg = 0;                                                             % probabilistic certificate found? (used in svm method) 
prob_flg_cvx = 0;                                                      % probabilistic back up certificate found? (used in cvx hull method) 

%% Essential arrays and initializations
x_cl = zeros(nx,simsteps+1,iC);
u_cl = zeros(nu,simsteps,iC);
cost_iter = zeros(1,iC); 
Xhat(1) = Xold;                                                         % initial estimate
Xhat_bck(1) = Xold;                                                  % backup set used in the cvx hull based method
x_cvx = [x_0'; x_ref']; 
lab_init = simul_con(x_cvx', Xnew);
x_net = x_cvx; 
y_net = lab_init'; 

%% Warm starting with some random data. 
% These can also be collected by applying random inputs to the system. For now, just using randomly sampled states.
%%%
% warm start data count
count_ws = 40;                                                         % good choices' tip:10 cvx, 40 for svm eps = .30,  60 for svm eps = .50 

% create some random points
rand_points = zeros(nx,count_ws);
for i = 1:count_ws
    rand_points(:,i) = -bold(1)*ones(nx,1) + 2*[bold(1)*rand(1); bold(3)*rand(1)]; 
end

% collect all the feasibility flags from the simulator
flags_ws = zeros(count_ws, size(Xnew.A,1));
lab_ws = zeros(1,count_ws); 
for i = 1:count_ws
    flags_ws(i,:) = simul_con(rand_points(:,i), Xnew);
    lab_ws(1,i)   = prod(flags_ws(i,:));  
end

%% Append the exploration data for a warm start
x_net = [x_net; rand_points']; y_net = [y_net; lab_ws'];
for i = 1:count_ws
        if lab_ws(1,i) == 1
            x_cvx = [x_cvx; rand_points(:,i)'];
        end
end

%% The main iteration loop starts here

for iter_count = 1:iC
 
    %%% getting Started with Closed Loop MPC Problem 
    x_cl(:,1,iter_count) = x_0;  
    options = sdpsettings('solver','gurobi','verbose',0);
    
    %% These are the constraint estimates 
    Chat_apnd = Xhat(iter_count).A; 
    bhat_apnd = Xhat(iter_count).b; 
    Chat = [Cold(1:4,:); Chat_apnd; Cold(5:6,:)]; 
    Dhat = [[0; 0; 0; 0]; zeros(size(Chat_apnd,1),1); [1; -1]]; 
    bhat = [bold(1:4); bhat_apnd; bold(5:6)]; 

    %% Checks needed for estimated constraints 
    % forming the terminal set with the estimated set
    [Xn,Pinf] = term_setRob(A,B,Chat,Dhat,bhat,Q,R,U,W,nx,nu);
    
    % check if it is empty
    term_flg = isEmptySet(Xn); 
    
    % feasibility check of the MPC problem from the starting state
    [matF, matG, matH, mat_c, dim_t] = obtain_matR(A,B,Chat,Dhat,bhat,Xn,N);
    [M, v, feas] = solveMatrices(A, B, N, W, dim_t, matF, mat_c, matH, matG, x_0, x_ref, Q, R, Pinf, 0); 
 
    if cvx_flag == 1 && term_flg == 0 && feas == 0 && iter_count > 1
        disp('****ROBUST CERTIFICATE OBTAINED****');
        rob_flg = 1;
        XhatcvxOut = Xhat(iter_count); 
        keyboard
    end
 
    count = 1;
    
    %% If not robust and paused, continue scaling estimates until checks are passed
    while (term_flg== 1 || feas == 1) 
        if cvx_flag == 0                                                                                   % inside SVM method
            disp('****SVM ESTIMATED SET BEING SCALED****');
            scale = (count/5)*step_size;                                                              % scaling strategy
            Xhat(iter_count) = intersect(Xold,scale*Xhat(iter_count));   

            % new estimates after scaling
            Chat_apnd = Xhat(iter_count).A; 
            bhat_apnd = Xhat(iter_count).b; 
            Chat = [Cold(1:4,:); Chat_apnd; Cold(5:6,:)]; 
            Dhat = [[0; 0; 0; 0]; zeros(size(Chat_apnd,1),1); [1; -1]]; 
            bhat = [bold(1:4); bhat_apnd; bold(5:6)]; 
    
            % forming terminal set again and checking emptiness
            [Xn,Pinf] = term_setRob(A,B,Chat,Dhat,bhat,Q,R,U,W,nx,nu);
            term_flg = isEmptySet(Xn); 
            
            % checking feasibility again with the scaled estimates
            [matF, matG, matH, mat_c, dim_t] = obtain_matR(A,B,Chat,Dhat,bhat,Xn,N);
            [M, v, feas] = solveMatrices(A, B, N, W, dim_t, matF, mat_c, matH, matG, x_0, x_ref, Q, R, Pinf, 0); 
            count = count + 1; 

        else                                                                                               % inside cvx hull method 
            
            if count == 1
                scale = 1;
            else
                disp('**** BACKUP SET BEING SCALED****');
                scale = (count/5)*step_size;                                                    % scaling strategy
            end
            
            % new backup and the constraints estimates
            Xhat_bck(iter_count) = intersect(Xold,scale*Xhat_bck(iter_count)); 
            Chat_apnd = Xhat_bck(iter_count).A; 
            bhat_apnd = Xhat_bck(iter_count).b; 
            Chat = [Cold(1:4,:); Chat_apnd; Cold(5:6,:)]; 
            Dhat = [[0; 0; 0; 0]; zeros(size(Chat_apnd,1),1); [1; -1]]; 
            bhat = [bold(1:4); bhat_apnd; bold(5:6)]; 
            
            % terminal set form and emptiness check
            [Xn,Pinf] = term_setRob(A,B,Chat,Dhat,bhat,Q,R,U,W,nx,nu);
            term_flg = isEmptySet(Xn); 
            
            % feasibility check with the new backup estimate
            [matF, matG, matH, mat_c, dim_t] = obtain_matR(A,B,Chat,Dhat,bhat,Xn,N);
            [M, v, feas] = solveMatrices(A, B, N, W, dim_t, matF, mat_c, matH, matG, x_0, x_ref, Q, R, Pinf, 0);
            count = count + 1; 
            
        end
    
    end

    %% Run an RMPC Along the Iteration (This will be feasible now after all scaling measures) 
    % solve MPC problem, updating closed-loop information

    for i = 1:simsteps
        i % printing the time step
        % solve optimization problem. feasible for sure!
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

    %% Iteration is over. Now need to update the estimated constraints
    % Need to form all cases here 
    [lab_iter] = simul_con(x_cl(:,:,iter_count), Xnew);
 
    if cvx_flag == 1                                                                                          % cvx hull method 
        disp('******INSIDE CVX METHOD******')
        if rob_flg == 0   
            disp('******NOT ROBUST YET******') 
            for j = 1:size(lab_iter,2)
                 if lab_iter(1,j) == 1                                                                       % satisfies constraints
                    x_cvx  = [x_cvx; x_cl(:,j,iter_count)']; 
                 end
            end
            Xhat(iter_count+1) = Polyhedron('V',x_cvx);    
            
            % KEEPING AN SVM BACKUP INSTEAD OF ARBITRARY SCALING %%  
            if prob_flg_cvx == 0   
                if prod(lab_iter) == 1                                                                     % only if violations are seen, update 
                    disp('******NO NEW VIOLATIONS******') 
                    if Xnew.contains(Xhat_bck(iter_count)) ==1
                       disp('******SVM FOUND AN INNER APPROXIMATION!!!******') 
                    end
                    Xhat_bck(iter_count+1) = Xhat_bck(iter_count); 
                    count_succ = count_succ + 1; 
                    if count_succ == Nlb_probVer
                        disp('PROBABILISTIC CERTIFICATE OBTAINED FOR CVX CASE');
                        prob_flg_cvx = 1; 
                    end
                else
                    disp('******NEW VIOLATIONS SEEN******') 
                    count_succ = 0; 
                    x_net = [x_net; x_cl(:,:,iter_count)'];                                           % using all data
                    y_net = [y_net; lab_iter']; 
                   %% Train SVM Boundary and See Performance 
                    SVMModel = fitcsvm(x_net, y_net, 'KernelFunction','rbf');
                    [innerpX] = visualizeBoundary(x_net, y_net', SVMModel);
                    polInner = Polyhedron('V',innerpX); 
                    Xhat_bck(iter_count+1)  = intersect(Xold,polInner);

                end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
            else
                disp('******SVM BACKUP SET CERTIFIED FOR CVX CASE******')
                Xhat_bck(iter_count+1)  = Xhat_bck(iter_count);                             % stop updates. probabilistic backup certificate obtained 
                keyboard
            end

        else
            disp('******CVX AND ROBUST******')
            Xhat(iter_count+1) = Xhat(iter_count);                                                % robust guarantees hold;  
            keyboard
        end
     
    else
        disp('******INSIDE SVM METHOD******')                                              % SVM method
        if prob_flg == 0   
            if prod(lab_iter) == 1                                                                          % only if violations are seen, update 
                 disp('******NO NEW VIOLATIONS******') 
                 if Xnew.contains(Xhat(iter_count)) ==1
                       disp('******SVM FOUND AN INNER APPROXIMATION!!!******') 
                 end
                 Xhat(iter_count+1) = Xhat(iter_count); 
                 count_succ = count_succ + 1; 
                 if count_succ == Nlb_probVer
                    disp('PROBABILISTIC CERTIFICATE OBTAINED');
                    prob_flg = 1; 
                    keyboard
                 end
            else
                disp('******NEW VIOLATIONS SEEN******') 
                count_succ = 0; 
                x_net = [x_net; x_cl(:,:,iter_count)'];                                                % using all data
                y_net = [y_net; lab_iter']; 
                %% Train SVM Boundary and See Performance 
                SVMModel = fitcsvm(x_net, y_net, 'KernelFunction','rbf');
                [innerpX] = visualizeBoundary(x_net, y_net', SVMModel);
                polInner = Polyhedron('V',innerpX); 
                Xhat(iter_count+1)  = intersect(Xold,polInner);
            end
        
        else
            disp('******SVM AND PROB CERTIFIED******')
            Xhat(iter_count+1)  = Xhat(iter_count);                                                % stop updates. Certificate obtained 
            XhatsvmOut = Xhat(iter_count);
            keyboard
        end
        
    end
    
    iter_count

end

%% Save
% Run this part once the code stops running
% Save all the data once estimated constraint sets converge
save('allData_XYZ.mat');
% Then run the compare_code.m file for all plots and tables.
