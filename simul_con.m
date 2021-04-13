%% Simulator Code (Checks Feasibility w.r.t all constraints and gives out a flag)
% Monimoy Bujarbaruah

function [lab] = simul_con(x_trajData, Xnew)

    flags = zeros(size(Xnew.A,1),size(x_trajData,2)); lab = zeros(1,size(x_trajData,2)); 
    tol_con = 1e-5;                                                     % constraint violation tolerance 

    for i = 1:size(x_trajData,2)
        flags(:,i) = ((Xnew.A)*x_trajData(:,i) <= Xnew.b + tol_con);
        lab(1,i)   = prod(flags(:,i));  
    end

     % Return all the violation labels 
     
end
