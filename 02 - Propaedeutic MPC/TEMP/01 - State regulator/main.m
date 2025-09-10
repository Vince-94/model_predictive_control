%% MAIN
close all
clear all
clc


%% Plant

T = 0.01;

A = [9.5310 195.4713; 0 -5.129];
B = [-7.6354; 8.0736];
C = [-1 1];
D = 0;
x0 = [1; 1];

% discretization
sys = ss(A, B, C, D);
sys_real = c2d(sys, T, 'zoh');
[Ad_sys, Bd_sys, Cd_sys, Dd_sys] = ssdata(sys_real);



%% Requirements

Ts = 0.01;

A_p = [9.5310 195.4713; 0 -5.129];
B_p = [-7.6354; 8.0736];
C_p = [-1 1];
D_p = 0;

% discretization
sys = ss(A_p, B_p, C_p, D_p);
sys_dp = c2d(sys, Ts, 'zoh');
[Ad_p, Bd_p, Cd_p, Dd_p] = ssdata(sys_dp);

% constraints
x_bound = [-inf inf
           -inf inf];
u_bound = [-10 10];
du_bound = [-inf inf];



%% MPC design
mpc_par.A = Ad_p;
mpc_par.B = Bd_p;
mpc_par.C = Cd_p;
mpc_par.n = size(A,1);
mpc_par.p = size(C,1);
mpc_par.q = size(B,2);
mpc_par.Ts = Ts;
mpc_par.Np = 50;                    % prediction horizon
mpc_par.Q = diag([10 10]);                    % output weight
mpc_par.R = diag([1]);                      % input weight
mpc_par.x_c = x_bound;              % state constrints
mpc_par.u_c = u_bound;              % input constrints
mpc_par.du_c = du_bound;            % input rate constrints

mpc_design = mpc_mtx_fun(mpc_par);




%% Simulation
t_sim = 1;

open('mpc_regulator')
tic;
sim('mpc_regulator')
toc;



%% Plot
% figure(1); hold on; grid on;
% plot(x.time, x.data(:,1), 'r');
% 
% figure(2); hold on; grid on;
% plot(x.time, x.data(:,2), 'r');