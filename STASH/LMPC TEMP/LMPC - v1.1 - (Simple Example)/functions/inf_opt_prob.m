function J_inf = inf_opt_prob(x, u, Q, R)

%{
Input:
    - X, U = stored states and inputs at j-th iteration
    - Q, R = weighting matrices
%}


k_inf = length(x);


for i = 1:k_inf

    % start from the end to the start
    j = k_inf-i+1;
    
    %% Iteration Cost: J = sum_inf(x'Qx + u'Ru)
    if i == 1
        Cost(j) = x(:,j)'*Q*x(:,j);
    else
        Cost(j) = Cost(j+1) + x(:,j)'*Q*x(:,j) + u(:,j)'*R*u(:,j);
    end
    
end

J_inf = Cost;