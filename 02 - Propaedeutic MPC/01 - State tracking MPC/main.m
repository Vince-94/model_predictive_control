%% MAIN
close all
clear all
clc


%% Requirements
T = 0.01;
Ts = 0.1;

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];
C = [1 0 0 0; 0 1 0 0];
D = zeros(size(C,1), size(B,2));

% discretization
sys = ss(A, B, C, D);
sys_dp = c2d(sys, Ts, 'zoh');
[Ad_p, Bd_p, Cd_p, Dd_p] = ssdata(sys_dp);

% reference
x0 = [0; 0; 0; 0];
ref = [10; 6; 0; 0];

% constraints
x_bound = [-inf inf
           -inf inf
           -inf inf
           -inf inf];
u_bound = [-inf inf
           -inf inf];
du_bound = [-inf inf
            -inf inf];



%% MPC design
mpc_param.A = Ad_p;
mpc_param.B = Bd_p;
mpc_param.C = Cd_p;
mpc_param.n = size(A,1);
mpc_param.p = size(C,1);
mpc_param.q = size(B,2);
mpc_param.Ts = Ts;
mpc_param.Np = 20;
mpc_param.Nc = 10;
mpc_param.Q = diag([1 1 1 1]);
mpc_param.R = diag([1 1]);
mpc_param.x_c = x_bound;
mpc_param.u_c = u_bound;
mpc_param.du_c = du_bound;

% mpc design
mpc = mpc_design(mpc_param);



%% Simulation
t_sim = 20;
open('mpc_tracking_state');

tic;
sim('mpc_tracking_state')
toc;


%% Plot
figure(1); hold on; grid on;
    title('Trajectory');
    plot(x.data(:,1), x.data(:,2), 'k');
    plot(x0(1), x0(2), 'ok');
    plot(ref(1), ref(2), 'xk');

figure(2);
subplot(3,2,1); hold on; grid on;
    title('x');
    plot(x.time, x.data(:,1), 'k');
subplot(3,2,2); hold on; grid on;
    title('y');
    plot(x.time, x.data(:,2), 'k');
subplot(3,2,3); hold on; grid on;
    title('v_x');
    plot(x.time, x.data(:,3), 'k');
subplot(3,2,4); hold on; grid on;
    title('v_y');
    plot(x.time, x.data(:,4), 'k');
subplot(3,2,5); hold on; grid on;
    title('a_x');
    plot(u.time, u.data(:,1), 'k');
subplot(3,2,6); hold on; grid on;
    title('a_y');
    plot(u.time, u.data(:,2), 'k');


