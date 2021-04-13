%% True Problem's Feasibility Check. If not, then no point
% Monimoy Bujarbaruah

%%
clear all
close all
clc
yalmip 'clear'

%% MPC Controller Parameters
[A,B,Cold, Dold, bold, Xold, Cnew, Dnew, bnew, Xnew, U, W, x_0, Q,R,N, x_ref,simsteps,nx, nu] = sys_load();

%% Find invariant set
[Xn,Pinf] = term_setRob(A,B,Cnew,Dnew,bnew,Q,R,U,W,nx,nu);
term_flg = isEmptySet(Xn); % must be 0 

keyboard

%% Forming matrices appearing in the optimization problem
% Checking with the new (full) constraints
[matF, matG, matH, mat_c, dim_t] = obtain_matR(A,B,Cnew,Dnew,bnew,Xn,N);

%% Check Feasibility
[M, v, feas] = solveMatrices(A, B, N, W, dim_t, matF, mat_c, matH, matG, x_0, x_ref, Q, R, Pinf, 0);
% feas must be 0