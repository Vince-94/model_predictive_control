function J_inf = cost_to_go(x, u, ref, Q, R)

%{
Input:
    - X, U = stored states and inputs at j-th iteration
    - Q, R = weighting matrices
%}


k_inf = length(x);


for i = 1:k_inf

    % start from the final state
    j = k_inf-i+1;
    
    
    %% Iteration Cost: J = sum_inf(x'Qx + u'Ru)
    if i == 1
        Cost(j) = (x(:,j)-ref)'*Q*(x(:,j)-ref);
    else
        Cost(j) = Cost(j+1) + (x(:,j)-ref)'*Q*(x(:,j)-ref) + u(:,j)'*R*u(:,j);
    end
    
end

J_inf = Cost;