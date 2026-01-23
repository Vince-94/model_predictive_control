function [X_LMPC, U_LMPC, X_stack, U_stack, J_inf] = lmpc(SS_0, U_0, Qfun_0, x0, xf, map, A, B, X_c, U_c, Np, Q, R);


%% Initialization
nx = length(A);


% Stacks initialization
X_stack{1} = SS_0;
U_stack{1} = U_0;
J_inf{1} = Qfun_0;



%% LMPC loop
% loop initialization
j = 1;              % initializing iteration
ExitFlag_1 = 0;      % initializing exit flag

while (ExitFlag_1 == 0)
    
    %% SS and Qfun polyhedron
    SSQfun = Polyhedron([SS_0', Qfun_0']);
%     SSQfun.computeHRep();
    SSQfun.computeVRep();
    SS = SSQfun.V(:,1:nx)';
    Qfun = SSQfun.V(:, end)';
    
    
    %% MPC loop
    % loop initialization
    X_LMPC = x0;                % initializing inital condition
    k = 1;                      % initializing sampling time
    ExitFlag_2 = 0;              % initializing exit flag
    
    
    while (ExitFlag_2 == 0)
        clc
        fprintf('Time step: k = %d, Iteration: j = %d, Initial Iteration Cost: J_iter(0) = %13.15f\n', [k, j, J_inf{j}(1)]);
        
        %% Finite Optimization Problem
        [u] = fin_opt_prob(X_LMPC(:,k), xf, Np, Q, R, Qfun, SS, A, B, X_c, U_c, map);
        
        %% Prediction
        U_LMPC(:,k) = u;
        X_LMPC(:,k+1) = A*X_LMPC(:,k) + B*U_LMPC(:,k);
        
        % plot
        figure(1); grid on; hold on;
            plot(X_LMPC(1,:), X_LMPC(2,:), 'r-o')       % trajectory
        
        
        %% Exit condition
        if (abs(X_LMPC(:,k+1) - xf) < 20e-2)
            ExitFlag_2 = 1;
        end
        
        
        %% Update sample
        k = k+1;
        
    end

    
    %% Store iteration data
    X_stack{j+1} = X_LMPC;
    U_stack{j+1} = U_LMPC;
    J_inf{j+1} = inf_opt_prob(X_LMPC, U_LMPC, Q, R);
    
    
    %% Update Safe Set and Q-function
    SS = [SS, X_LMPC];
    Qfun = [Qfun, J_inf{j+1}];
    
    
    %% Exit condition
    if (abs(J_inf{j}(1) - J_inf{j+1}(1)) < 10e-2)
        ExitFlag_1 = 1;
    end
    
    
    %% Update iteration
    j = j + 1;
    
end

end