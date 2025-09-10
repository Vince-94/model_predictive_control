%% MAIN
close all
clear all
clc


%% Plant
T = 0.001;

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];
C = [1 0 0 0; 0 1 0 0];
D = zeros(size(C,1), size(B,2));
x0 = [1; 1; 0; 0];
y0 = C*x0;

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
ref = [7; 1];

% constraints
x_bound = [-Inf Inf
           -Inf Inf
           -2 2
           -2 2];
u_bound = [-1 1
           -1 1];
du_bound = [-1 1
            -1 1];
        
% obstacles
sx = [0; 2; 4; 6; 8; 10; 12; 14; 16; 20];
sy = [0; 3; 7; 10];



%% Map parametrization
obs_param = obs_map_fun(sx, sy);
obs_param.r_act = 2;
obs_param.r_safe = 0.5;



%% Waypoints
wp.points = [[1;8] [5;8] [5;2] [9;2] [9;8] [13;8] [13;2] [18;2] [18;9]];
% wp.points = [[3;8] [7;2] [11;8] [15;2] [18;9]];
% wp.points = [[1;8] [5;8] [6;2] [9;2]];
wp.tol = 0.5;



%% MPC design
mpc_par.A = Ad_p;
mpc_par.B = Bd_p;
mpc_par.C = Cd_p;
mpc_par.n = size(Ad_p,1);
mpc_par.p = size(Cd_p,1);
mpc_par.q = size(Bd_p,2);
mpc_par.Ts = Ts;
mpc_par.Np = 30;                    % prediction horizon
mpc_par.Nc = 30;                     % control horizon
mpc_par.Q = diag([1 1]);                    % output weight
mpc_par.R = diag([1 1]);                      % input weight
mpc_par.x_c = x_bound;              % state constrints
mpc_par.u_c = u_bound;              % input constrints
mpc_par.du_c = du_bound;            % input rate constrints

mpc_design = mpc_mtx_fun(mpc_par);



%% Simulation
t_sim = 200;
open('mpc_v2');

tic;
sim('mpc_v2')
toc;



%% Plot
figure(1); hold on;
    plot(y.data(:,1), y.data(:,2));
    plot(x0(1), x0(2), 'ro');
    
    for i=1:length(wp.points)
        plot(wp.points(1,i), wp.points(2,i), 'rx')
    end

    for i=1:length(obs_param.Sx)
        plot([obs_param.Sx(i,1);obs_param.Sx(i,1)], [obs_param.Sx(i,2);obs_param.Sx(i,3)], 'k')
    end
    
    for i=1:length(obs_param.Sy)
        plot([obs_param.Sy(i,2);obs_param.Sy(i,3)], [obs_param.Sy(i,1);obs_param.Sy(i,1)], 'k')
    end 

    
%% Save plot
name = sprintf('T=%.3f, Ts=%.2f, Np=%i, Nc=%i, Q=[%i,%i], R[%i,%i], u_bound=[%i,%i].jpeg', T, Ts, mpc_par.Np, mpc_par.Nc, mpc_par.Q(1,1), mpc_par.Q(2,2), mpc_par.R(1,1), mpc_par.R(2,2), u_bound(1,1), u_bound(1,2))
saveas(figure(1), name)


