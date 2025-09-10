%% MAIN
close all
clear all
clc

addpath('functions','-end');


%% Requirements

T = 0.01;
x0 = [1; 1; 0; 0];
ref = [18; 9; 0; 0];

Ts = 0.1;
Np = 30;
Nc = 30;
Q = diag([1 1 1 1]);
R = diag([1 1]);
[pred_model] = point_lin_kin(Ts);

r_act = 2;
r_safe = 0.25;
r_tol = 0.2;



%% MPC design
mpc_param.Ts = Ts;
mpc_param.Np = Np;
mpc_param.Nc = Nc;
mpc_param.Q = Q;
mpc_param.R = R;
mpc_param.A = pred_model.Ad;
mpc_param.B = pred_model.Bd;
mpc_param.C = pred_model.Cd;
mpc_param.n = size(pred_model.Ad,1);
mpc_param.p = size(pred_model.Cd,1);
mpc_param.q = size(pred_model.Bd,2);
mpc_param.x_c = pred_model.x_bound;
mpc_param.u_c = pred_model.u_bound;
mpc_param.du_c = pred_model.du_bound;

[mpc] = mpc_design(mpc_param);



%% Obstacle avoidance

% Map definition
[map] = obs_design(r_act, r_safe);

% Waypoints definition
[wp] = waypoint_design(r_tol);



%% Simulation
t_sim = 120;
open('mpc_v0_1');

tic;
sim('mpc_v0_1')
toc;



%% Plot
figure(1); hold on; grid on;
    title('Trajectory');
    plot(x0(1), x0(2), 'ok');
    plot(ref(1), ref(2), 'xk');
    plot(x.data(:,1), x.data(:,2), 'k');
    % plot environment
%     for i=1:length(map.Sx)
%         plot([map.Sx(i,1);map.Sx(i,1)], [map.Sx(i,2);map.Sx(i,3)], 'k')
%     end
%     for i=1:length(map.Sy)
%         plot([map.Sy(i,2);map.Sy(i,3)], [map.Sy(i,1);map.Sy(i,1)], 'k')
%     end 
    for i=1:length(wp.points)-1
        plot(wp.points(1,i), wp.points(2,i), 'ro')
    end
    
    
figure(2);
subplot(3,2,1); hold on; grid on;
    title('x');
    plot(x.time, x.data(:,1), 'k');
    plot(target.time, target.data(:,1), 'r');
subplot(3,2,2); hold on; grid on;
    title('y');
    plot(x.time, x.data(:,2), 'k');
    plot(target.time, target.data(:,2), 'r');
subplot(3,2,3); hold on; grid on;
    title('v_x');
    plot(x.time, x.data(:,3), 'k');
    plot(target.time, target.data(:,3), 'r');
subplot(3,2,4); hold on; grid on;
    title('v_y');
    plot(x.time, x.data(:,4), 'k');
    plot(target.time, target.data(:,4), 'r');
subplot(3,2,5); hold on; grid on;
    title('a_x');
    plot(u.time, u.data(:,1), 'k');
subplot(3,2,6); hold on; grid on;
    title('a_y');
    plot(u.time, u.data(:,2), 'k');


%% Save data
% save('feasible_solution.mat', 'x', 'u')
save('reference.mat', 'trajectory')