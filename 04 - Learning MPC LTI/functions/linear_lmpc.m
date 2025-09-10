function [X, U, J_inf, SS] = linear_lmpc(X0, U0, x0, A, B, X_c, U_c, Np, Q, R)


%% Stack initialization
X{1} = X0;
U{1} = U0;
J_inf{1} = inf_opt_prob(X0, U0, Q, R);


%% Polyedron initialization
SSQfun = Polyhedron([X{1}', J_inf{1}']);
SSQfun.computeHRep();
SSQfun.computeVRep();
SS = SSQfun.V(:,1:length(A))';
Q_fun = SSQfun.V(:, end)';



%% LMPC loop
j = 1;               % initializing iteration
ExitFlag_1 = 0;      % initializing exit flag

while (ExitFlag_1 == 0)
    
    X_LMPC = x0;                % initializing inital condition
    k = 1;                      % initializing sampling time
    ExitFlag_2 = 0;              % initializing exit flag
    
    %% MPC loop
    while (ExitFlag_2 == 0)
        clc
        fprintf('Time step: k = %d, Iteration: j = %d, Initial Iteration Cost: J_iter(0) = %13.15f\n', [k, j, J_inf{j}(1)]);
        
        
        %% Finite Horizon Optimization Problem
        U_LMPC(:,k) = linear_mpc(X_LMPC(:,k), Np, Q, R, Q_fun, SS, A, B, X_c, U_c);
%         U_LMPC(:,k) = U_LMPC(:,k) + 0.01*rand;
        
        %% Prediction
        X_LMPC(:,k+1) = A*X_LMPC(:,k) + B*U_LMPC(:,k);
%         X_LMPC(:,k+1) = X_LMPC(:,k+1) + 0.001*[rand; rand];
        
         
        %% Exit condition
        if (abs(X_LMPC(:,k+1)) <= 10e-4)
            ExitFlag_2 = 1;
        end
        
        %% Update sample
        k = k+1;
        
    end
    
    
    %% Store iteration data
    X{j+1} = X_LMPC;
    U{j+1} = U_LMPC;
    J_inf{j+1} = inf_opt_prob(X_LMPC, U_LMPC, Q, R);
    
    
    %% Update Safe Set and Q-function
    SS = [SS, X_LMPC];
    Q_fun = [Q_fun, J_inf{j+1}];
    
    
    %% Exit condition
    if (abs(J_inf{j}(1) - J_inf{j+1}(1)) < 10e-4)
        ExitFlag_1 = 1;
    end
    
    
    %% Update iteration
    j = j + 1;
    
end


end