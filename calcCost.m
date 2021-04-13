%% Calculate the MPC Closed Loop Cost
% Monimoy Bujarbaruah and Charlott Vallon

function cost = calcCost(x, u, Q, R, Pinf, x_ref)

    cost = 0;
    for i = 1:length(u)
        cost = cost + u(i)*R*u(i) + (x(:,i)-x_ref)'*Q*(x(:,i)-x_ref);
    end
    cost = cost + (x(:,end)-x_ref)'*Pinf*(x(:,end)-x_ref);

end