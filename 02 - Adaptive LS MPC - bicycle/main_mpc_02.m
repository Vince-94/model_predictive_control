%% MAIN
close all
clear all
clc


addpath('functions','-end');
addpath('data','-end');


%% Reference
load('reference');

ref = trajectory;




%% =================================================================
%  ============================ Control ============================
%  =================================================================

T = 0.01;

T_ctrl = 0.5;
x0 = [1; 1; pi/2];
x_ref = [18; 8.7; 0];
v_des = 1;

% Constraints
x_min = [-inf; -inf; -inf];
x_max = [inf; inf; inf];
u_min = [-0.5];
u_max = [0.5];
du_min = [-0.5];
du_max = [0.5];


% MPC control
mpc_ctrl_param.Ts = T_ctrl;
mpc_ctrl_param.nx = 3;
mpc_ctrl_param.nu = 1;
mpc_ctrl_param.ny = 2;
mpc_ctrl_param.Np = 50;
mpc_ctrl_param.Nc = 5;
mpc_ctrl_param.Q = diag([10 10]);
mpc_ctrl_param.R = diag([1]);
mpc_ctrl_param.x_c = [x_min, x_max];
mpc_ctrl_param.u_c = [u_min, u_max];
mpc_ctrl_param.du_c = [du_min, du_max];




%% ==================================================================
%  =========================== Simulation ===========================
%  ==================================================================
t_sim = Inf;

open('mpc_02_sim');

tic;
sim('mpc_02_sim')
toc;



%% Plot
plot_trajectory(ref, 'g');
plot_trajectory(x, 'k');

figure(1); grid on; hold on;
plot(target.data(end,1), target.data(end,2), 'xr');


return