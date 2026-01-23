%% MAIN
clc
clear all
close all

addpath('data','-end')
addpath('functions','-end')


%% Documentation
%{
    States: x = [x, y, vx, vy]
    Inputs: u = [ax, ay]

%}



%% System model
Ts = 0.1;
[A, B, X_c, U_c] = point_mass_kin(Ts);



%% Obstacle avoidance
r_act = 2;
r_safe = 0.25;
[map] = obs_design(r_act, r_safe);



%% MPC parameters
Np = 4;                        % prediction horizon
Q = diag([1 1 1 1]);            % state cost matrix
R = diag([1 1]);                % input cost matrix




%% ==================================================================
%  ============================== LMPC ==============================
%  ==================================================================
%{
    X0, U0 = initial feasible trajectory/input
    J_inf0 = initial iteration cost (Q_function)

    X_stack{j}, U_stack{j} = trajectory/input stack
    J_inf{j} = iteration cost stack (Q_function)

    X_LMPC, U_LMPC = 
    
%}


%% Load the feasible solution
load('data/feasible_solution.mat')
n_points = 100;

% state
X0 = x.data';
space = round(length(X0)/n_points);
X0 = X0(:, 1:space:end);

% input
U0 = u.data';
space = round(length(U0)/n_points);
U0 = U0(:, 1:space:end);

% initial and final state
x0 = [X0(1,1); X0(1,2); 0; 0];                      % initial state
xf = [X0(1,end); X0(2,end); 0; 0];                  % final state

% initial infinite cost
J_inf0 = inf_opt_prob(X0, U0, Q, R);



%% Plot
figure(1); grid on; hold on;
    title('Trajectories');
    xlabel('$$x_1$$','interpreter','latex','fontsize',15);
    ylabel('$$x_2$$','interpreter','latex','fontsize',15);
    plot(X0(1,:), X0(2,:), 'k-o');      % feasible solution
    plot(x0(1), x0(2), '^');            % starting point
    plot(xf(1), xf(2), 's');            % final goal




%% LMPC
[X, U, X_stack, U_stack, J_inf] = lmpc(X0, U0, J_inf0, x0, xf, map, A, B, X_c, U_c, Np, Q, R);



%% Perfermance
iters = [1:length(X_stack)];

for i = 1:length(X_stack)
    J_iter0(i) = J_inf{i}(1);
    fprintf('Iteration cost: J_iter0(%d) = %13.15f\n', [i, J_iter0(i)]);
end

figure(2); hold on; grid on;
    plot(iters-1, J_iter0(:), 'r-o');       % Iteration cost at starting trajectory


