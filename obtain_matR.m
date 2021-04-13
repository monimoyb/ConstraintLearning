%% Forming Matrices for Robust Optimization (See https://www.sciencedirect.com/science/article/pii/S0005109806000021)

function [matF, matG, matH, mat_c, dim_t] = obtain_matR(A,B,C,D,b,Xn,N)
    
    dim_t = size(C,1)*N + size(Xn.A,1);
    nx = size(A,2);
    nu = size(B,2);
    capA = eye(nx);

    for k=1:N
        capA = [capA;A^k];
    end

    matE(:,:,1) = [eye(nx),zeros(nx,nx*(N-1))];
    capE = [zeros(nx,nx*N);matE(:,:,1)];

    for k=2:N
        matE(:,:,k) = [A^(k-1),matE(:,1:(end-nx),k-1)];
        capE = [capE;matE(:,:,k)];
    end

    capB = capE*kron(eye(N),B);
    capC = blkdiag(kron(eye(N),C),Xn.A);
    capD = [kron(eye(N),D);zeros(dim_t-size(kron(eye(N),D),1),nu*N)];

    %% Batch matrices for simulation check
    Aw_batch = capE(nx+1:end,:);
    Bu_batch = capB(nx+1:end,:);
    A_batch  = capA(nx+1:end,:);

    matF = capC*capB + capD;
    matG = capC*capE;
    matH = -capC*capA;
    mat_c = [kron(ones(N,1),b); Xn.b];

end