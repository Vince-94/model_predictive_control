%% MAIN
close all
clear all
clc

addpath('mpc_functions','-end')
addpath('curvilinear_functions','-end')
addpath('utils','-end')
addpath('plant','-end')


%% =================================================================
%  ============================ Plant ============================
%  =================================================================

% UGV = [vx, vy, wz, X, Y, psi];
x0 = [1; 0; 0; 0; 0; 0];
T = 0.01;

figure(1); hold on; grid on;
title('Trajectory');
plot(x0(4), x0(5), 'or');
arrow(x0(4), x0(5), x0(6), 'r');

% vehicle data
sys_model.lf = 0.75;
sys_model.lr = 0.75;
sys_model.L = sys_model.lf+sys_model.lr;
sys_model.W = 0.8;
sys_model.m = 800;
sys_model.Iz = sys_model.m*(sys_model.L^2 + sys_model.W^2)/12;
sys_model.Cf = 10000;
sys_model.Cr = 10000;




%% =================================================================
%  ============================ Track ==============================
%  =================================================================

track_n = 5;
W_track = 2;
[PointAndTangent, x0_track] = track_fun(track_n, W_track);




%% ==================================================================
%  ======================= MPC initialization =======================
%  ==================================================================

T_ctrl = 0.1;
method = 0;

% initial condition
X_plant1 = x0;
% X_plant2 = x0;
X_plant2 = [x0(4); x0(5)];
[s0, ey0, epsi0] = global2local_fun(x0(4:6), x0_track, PointAndTangent, W_track);
X_local2 = [x0(1); x0(2); x0(3); s0; ey0; epsi0];


% target
target = [PointAndTangent(end,1); PointAndTangent(end,2); PointAndTangent(end,3)];
[s_target, ey_target, epsi_target] = global2local_fun(target, x0(4:6), PointAndTangent, W_track);
ref = [0; 0; 0; s_target-s0; ey_target; epsi_target];




%% ==================================================================
%  ============================== MPC ===============================
%  ==================================================================


% MPC tuning
Np = 10;
Nc = Np;
Q = diag([1 10 10 1 100 100]);
R = diag([10 10]);
S = diag([0 0]);

k = 1;
flag = 0;
t_sim = 10;

U_old = zeros(Nc*2,1);

while (flag == 0)
    
    clc
    t = T_ctrl*k;
    
    %% Print info
    fprintf('Sample: k = %d, time: t = %.2f\n', [k, t]);
    
    
    %% MPC
    U_MPC(:,k) = [0.35; 0];
    fprintf('Input =            [delta_f = %.1f°,  a = %.2fm/s^2]\n', [U_MPC(1,k)*180/pi, U_MPC(2,k)]);


    %% Plant
    x1 = X_plant1(:,k);
    x2 = X_local2(:,k);
    
    u = U_MPC(:,k);
    
    for i=1:T_ctrl/T
        % plant 1
        [x1_plus, alpha] = bicycle_plant_discrete(x1, u, sys_model, T);
        x1 = x1_plus;
        
        % plant 2
        model2 = bicycle_pred_model(x2, u, T, PointAndTangent);
        x2_plus = model2.Ad*x2 + model2.Bd*u + model2.Kd;
        x2 = x2_plus;
        
    end
    
    eigenvalues = eig(model2.Ad)';
    
    X_plant1(:,k+1) = x1;
    X_local2(:,k+1) = x2;
    
    [x_plant2, y_plant2] = local2global_fun(X_local2(4:6,k+1), [s0; ey0; epsi0], PointAndTangent);
    X_plant2(:,k+1) = [x_plant2; y_plant2];

    
    %% Update sample
    k = k+1;
    u_old = U_MPC(:,k-1);
    
    figure(1); hold on; grid on; axis equal;
    
    p2 = plot(X_plant1(4,k), X_plant1(5,k), '.r', 'LineWidth', 0.5);
%     arrow(X_plant1(4,k), X_plant1(5,k), X_plant1(6,k), 'r');
    
%     p3 = plot(X_plant2(4,k), X_plant2(5,k), 'ob', 'LineWidth', 0.5);
%     arrow(X_plant2(4,k), X_plant2(5,k), X_plant2(6,k), 'b');
    p3 = plot(X_plant2(1,k), X_plant2(2,k), '.b', 'LineWidth', 0.5);

    

    %% Exit condition
    if (t == t_sim)
        fprintf('Simulation time reached \n');
        flag = 1;
        t_fin = T_ctrl*(k-1);
    end
    
end

