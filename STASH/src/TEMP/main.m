%% MAIN
close all
clear all
clc


%% Plant
T = 0.01;

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];
C = [1 0 0 0; 0 1 0 0];
D = zeros(size(C,1), size(B,2));
x0 = [2; 1; 0; 0];

% discretization
sys = ss(A, B, C, D);
sys_real = c2d(sys, T, 'zoh');
[Ad_sys, Bd_sys, Cd_sys, Dd_sys] = ssdata(sys_real);



%% Requirements
Ts = 0.1;

A_p = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B_p = [0 0; 0 0; 1 0; 0 1];
C_p = [1 0 0 0; 0 1 0 0];
D_p = zeros(size(C_p,1), size(B_p,2));

% discretization
sys = ss(A_p, B_p, C_p, D_p);
sys_dp = c2d(sys, Ts, 'zoh');
[Ad_p, Bd_p, Cd_p, Dd_p] = ssdata(sys_dp);

% reference
ref = [8; 10];

% constraints
x_bound = [-inf inf
           -inf inf
           -inf inf
           -inf inf];
u_bound = [-5 5
           -5 5];
du_bound = [-inf inf
            -inf inf];


%% MPC design
mpc_par.A = Ad_p;
mpc_par.B = Bd_p;
mpc_par.C = Cd_p;
mpc_par.n = size(A,1);
mpc_par.p = size(C,1);
mpc_par.q = size(B,2);
mpc_par.Ts = Ts;
mpc_par.Np = 20;
mpc_par.Nc = 10;
mpc_par.Q = diag([1 1]);
mpc_par.R = diag([1 1]);
mpc_par.x_c = x_bound;
mpc_par.u_c = u_bound;
mpc_par.du_c = du_bound;

mpc_design = mpc_mtx_fun(mpc_par);



%% Simulation
t_sim = 20;
open('mpc_control_horizon');

tic;
sim('mpc_control_horizon')
toc;


%% Plot
figure(1); hold on;
    plot(y.data(:,1), y.data(:,2));
    plot(x0(1), x0(2), 'ok');
    plot(ref(1), ref(2), 'xr');
