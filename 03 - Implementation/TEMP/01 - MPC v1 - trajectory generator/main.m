%% MAIN
close all
clear all
clc


%% Requirements
Ts = 0.1;

% system
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];
C = [1 0 0 0; 0 1 0 0];
% C = eye(4);
D = zeros(size(C,1), size(B,2));

% initial conditions
x0 = [2; 1; 0; 0];

% reference
target = [18; 6];
vx_des = 0;
vy_des = 0;
ref = [target];

% design Constraints
x_bound = [-100 100];
u_bound = [-100 100];
du_bound = [-50 50];



%% MPC design

n = size(A,1);
q = size(B,2);
p = size(C,1);

% prediction model
sys = ss(A, B, C, D);
sys_d = c2d(sys, Ts, 'zoh');
[Ad, Bd, Cd, Dd] = ssdata(sys_d);

% horizon
Np = 50;

% weighing matrices
Q = diag([10 10]);
R = diag([1 1]);



%% MPC functions

% matrices initialization
U0 = zeros(q*Np, 1);
A_bar = zeros(n*(Np+1), n);
B_bar = zeros(n*(Np+1), q*Np);
C_bar = zeros(p*(Np+1), n);
D_bar = zeros(p*(Np+1), q*Np);


% prediciton matrices
for i = 1:Np+1
    for j = 1:Np
        
        A_bar((i-1)*n+1:i*n, 1:n) = Ad^(i-1);
        C_bar((i-1)*p+1:i*p, 1:n) = Cd*Ad^(i-1);
        
        if j < i
            B_bar((i-1)*n+1:i*n, (j-1)*q+1:j*q) = Ad^(i-j-1)*Bd;
            D_bar((i-1)*p+1:i*p, (j-1)*q+1:j*q) = Cd*Ad^(i-j-1)*Bd;
        else
            B_bar(i, j) = 0;
            D_bar(i, j) = 0;
        end
    end
end


% objective function matrices
Q_cell = repmat({Q}, Np+1, 1);
Q_bar = blkdiag(Q_cell{:});
R_cell = repmat({R}, Np, 1);
R_bar = blkdiag(R_cell{:});


% constraint matrices
X_bound = ones(n*(Np+1), 1)*x_bound;
U_bound = ones(q*Np, 1)*u_bound;
dU_bound = ones(q*Np, 1)*du_bound;





%% Simulation

t_sim = 10;
open('mpc_two_layers')
tic;
sim('mpc_two_layers')
toc;


figure(1); hold on; grid on;
plot(y.data(:,1), y.data(:,2));
plot(x0(1), x0(2), 'ok');
plot(target(1), target(2), 'xr');
