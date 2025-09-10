%% MAIN
close all
clear all
clc


%% Requirements
Ts = 0.01;

% system
A = [9.5310 195.4713; 0 -5.129];
B = [-7.6354; 8.0736];
C = [-1 1];
D = 0;

% initial conditions
x0 = [0; 0];

% reference
ref = 2;

% design Constraints
x_bound = [-inf inf;
           -inf inf];
u_bound = [-inf inf];
du_bound = [-inf inf];



%% MPC design

n = size(A,1);
q = size(B,2);
p = size(C,1);

% prediction model
sys = ss(A, B, C, D);
sys_d = c2d(sys, Ts, 'zoh');
[Ad, Bd, Cd, Dd] = ssdata(sys_d);

% horizon
Np = 20;

% weighing matrices
Q = 10;
R = 1;



%% MPC function

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
X_bound = repmat(x_bound, Np+1, 1);
U_bound = repmat(u_bound, Np, 1);
dU_bound = repmat(du_bound, Np, 1);



%% Simulation

t_sim = 2;
open('mpc_tracking');

tic;
sim('mpc_tracking')
toc;



%% Plot
figure(1); hold on; grid on;
plot(y.time, y.data, 'r');
plot(y_pred.time, y_pred.data, 'k');
