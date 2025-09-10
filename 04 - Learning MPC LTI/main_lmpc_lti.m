%% MAIN
clc
clear all
close all

addpath('functions','-end')


%% System
A = [1 1; 0 1];
B = [0; 1];

% dimension
nx = size(A,1);
nu = size(B,2);

% initial state
x0 = [-3.95; -0.05];


% state constraint polihedron
x_min = [-5; -5];
x_max = [5; 5];
A_x = [-1*eye(nx);
       1*eye(nx)];
b_x = [-x_min; x_max];
X_c = Polyhedron(A_x, b_x);

% input constraint polihedron
u_min = -1;
u_max = 1;
A_u = [-1*eye(nu);
       1*eye(nu)];
b_u = [-u_min; u_max];
U_c = Polyhedron(A_u, b_u);



%% LQR
Q = diag([1 1]);
R = 1;

flag = 10;
k = 1;
x(:,1) = x0;
ExitFlag = 0;

figure(1); hold on;
plot(x(1,1), x(2,1), 'ok');

while ExitFlag == 0
    [K,P] = dlqr(A, B, Q, R);
    u(k) = -K*x(:,k);
    
    k = k+1;
    
    x(:,k) = A*x(:,k-1) + B*u(k-1);
    
    figure(1); hold on;
    plot(x(1,:), x(2,:), 'o-k');
    
    figure(2); hold on;
    plot(k, x(1,k), 'ob');
    plot(k, x(2,k), 'or');
    
    if x(:,k) <= 10e-4
        ExitFlag = 1;
    end
end
    
close all;




    
%% ==================================================================
%  ============================== LMPC ==============================
%  ==================================================================

%% Feasible solution
X0 = x;
U0 = u;


%% MPC
Ts = 0.1;
Np = 4;
Q = diag([1 1]);
R = 1;


%% LMPC Initialization

% Stack initialization
X{1} = X0;
U{1} = U0;
J_inf{1} = inf_opt_prob(X0, U0, Q, R);

% Polyedron initialization
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





%% ===================================================================
%  =========================== Performance ===========================
%  ===================================================================

SS = Polyhedron(SS');

for i = 1:length(X)-1
    
    % Trajectories
    figure(1); hold on; grid on;
    title('Trajectories');
    xlabel('$$x_1$$','interpreter','latex','fontsize',15);
    ylabel('$$x_2$$','interpreter','latex','fontsize',15);
    plot(SS, 'color', 'lightyellow');
    plot(X{i}(1,:), X{i}(2,:), '-o', 'LineWidth', 0.75);
    
    % Performance
    figure(2); hold on; grid on;
    title('Performance');
    xlabel('$$k$$','interpreter','latex','fontsize',15);
    ylabel('$$J_{\infty}$$','interpreter','latex','fontsize',15);
    plot(find(J_inf{i}(:))-1, J_inf{i}(:), '-o');

end

