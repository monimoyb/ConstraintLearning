%% Defining System and Constraint Matrices 
% Monimoy Bujarbaruah

function [A,B,Cold,Dold,bold,Xold,Cnew,Dnew,bnew,Xnew,U,W, x_0, Q,R,N, x_ref,simsteps, nx, nu] = sys_load()

    %% Considering two states and one scalar input 
    A = [1, 1; 0, 1]; 
    B = [0;1];   
    nx = size(A,2); nu = size(B,2); 
    %% Weights, horizon and task duration
    Q =  10*eye(nx);
    R =   2*eye(nu);
    N = 4;
    simsteps = 10;                                                                                
    %% Constraints of the form -a<=x(i)<=a and -ulb<=u<=uub
    % Expressing in Cx+Du <=b format 
    % These are the known constraints
    Cold = [1 0; -1 0; 0 1; 0 -1; 0 0; 0 0]; 
    Dold = [0; 0; 0; 0; 1; -1]; 
    bold = [20;20;20;20;30;30]; 
    Xold = Polyhedron('A',Cold(1:4,:),'b',bold(1:4,:));
    U = Polyhedron('A',Dold(5:6,:),'b',bold(5:6,:));        
    %% Adding the unknown constraints here 
    % Adding 2 new hyperplanes 
    Cnew = [Cold(1:4,:); [1 1; 1 -1]; Cold(5:6,:)]; 
    Dnew = [0; 0; 0; 0; 0; 0; 1; -1]; 
    bnew = [bold(1:4); [5; 5]; bold(5:6)];                     
    Xnew = Polyhedron('A',Cnew(1:6,:),'b',bnew(1:6,:));
    %% Defining Noise  Bounds 
    wub_true =  0.5;                                   % Upper bound
    wlb_true = -0.5;                                   % Lower bound 
    W = Polyhedron('lb',wlb_true*ones(nx,1),'ub',wub_true*ones(nx,1));
    %% Starting condition and reference
    x_0 = [-15; 15];    
    x_ref = [5; 0];    

end
