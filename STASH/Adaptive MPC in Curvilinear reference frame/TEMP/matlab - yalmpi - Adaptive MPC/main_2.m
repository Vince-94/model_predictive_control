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
sys_model.Cf = 2000;
sys_model.Cr = 2000;



%% =================================================================
%  ============================ Track ==============================
%  =================================================================

track_n = 2;
W_track = 2;
[PointAndTangent, x0_track] = track_fun(track_n, W_track);




%% ==================================================================
%  ======================= MPC initialization =======================
%  ==================================================================

T_ctrl = 0.1;

% initial condition
X = x0;
[s0, ey0, epsi0] = global2local_fun(x0(4:6), x0_track, PointAndTangent, W_track);
X_MPC = [x0(1); x0(2); x0(3); s0; ey0; epsi0];
u0 = [0; 0];

% target
target = [PointAndTangent(end,1); PointAndTangent(end,2); PointAndTangent(end,3)];
[s_target, ey_target, epsi_target] = global2local_fun(target, x0(4:6), PointAndTangent, W_track);
ref = [0; 0; 0; s_target-s0; ey_target; epsi_target];


% state constraints [vx, vy, wz, s, ey, epsi]
x_c = [0,           +3;
       -0.2,        +0.2;
       -50*pi/180,  +50*pi/180;
       0,           +s_target+0.5;
       -W_track,    +W_track;
       -10*pi/180,  +10*pi/180];


% input constraints [delta_f, a]
u_c = [-20*pi/180,  +20*pi/180;
       -0.5,        +0.5];



%% ==================================================================
%  ============================== MPC ===============================
%  ==================================================================

% MPC tuning
Np = 20;
Q = diag([1 10 10 1 10 10]);
R = diag([1 1]);

% loop parameters
t_sim = 25;         % maximum simulation time
k = 1;              % initial sample (k = 0)
flag = 0;           % initialize exit flag


while (flag == 0)
    
    t = T_ctrl*(k-1);
    %% Print Info
    clc
    fprintf('Sample: k = %d, time: t = %.2f\n', [k-1, t]);
    fprintf('Velocity =         [vx = %.2fm/s,  vy = %.2fm/s,   wz =   %.1f°/s] \n', [X_MPC(1,k), X_MPC(2,k), X_MPC(3,k)*180/pi]);
    fprintf('Curvilinear pose = [s =  %.2fm,    ey = %.2fm,     epsi = %.1f°] \n',   [X_MPC(4,k), X_MPC(5,k), X_MPC(6,k)*180/pi]);
    fprintf('Input =            [delta_f = %.1f°,  a = %.2fm/s^2]\n', [u0(1)*180/pi, u0(2)]);


    
    %% Prediction
    prediction_model = bicycle_pred_model(X_MPC(:,k), u0, T_ctrl, PointAndTangent);

    

    %% Optimization
    U_MPC(:,k) = adaptive_mpc(X_MPC(:,k), ref, Np, Q, R, prediction_model, x_c, u_c);
    if isnan(U_MPC(:,k))
        fprintf('Problem infeasible at k = %d \n', k);
        t_fin = t;
        U_MPC(:,k) = [];
        break;
    end
    

    %% Plant discrete
    x = X(:,k);
    figure(1);
    for i=1:T_ctrl/T
        [x(:,i+1), alpha] = bicycle_plant_discrete(x(:,i), U_MPC(:,k), sys_model, T);
        p1(i) = plot(x(4,i+1), x(5,i+1), '.k');
    end
    
    X(:,k+1) = x(:,end);
    alpha_f(k) = alpha;
    
    p2 = plot(X(4,k+1), X(5,k+1), 'or', 'LineWidth', 0.5);
    arrow(X(4,k+1), X(5,k+1), X(6,k+1), 'b');
    
    
    %% Conversion from global to curvilinear coordinates
    [s, ey, epsi, CompletedFlag] = global2local_fun(X(4:6,k+1), x0_track, PointAndTangent, W_track);
    X_MPC(:,k+1) = [X(1,k+1); X(2,k+1); X(3,k+1); s; ey; epsi];
    
    if CompletedFlag == 0
        flag = 1;
        t_fin = t;
        break;
    end
    

    
    %% Update sample
    k = k+1;
    u0 = U_MPC(:,k-1);
    
    
    
    %% Exit condition
    if (abs(X_MPC(4,k) - s_target) <= 10e-3 && X_MPC(5,k) <= 0.5)
        fprintf('Terget reached \n');
        t_fin = T_ctrl*(k-1);
        flag = 1;
    elseif (t == t_sim)
        fprintf('Simulation time reached \n');
        t_fin = T_ctrl*(k-1);
        flag = 1;
    end
    
    
end



%% ==================================================================
%  =========================== Plot =================================
%  ==================================================================

t_vec = linspace(0, t_fin, (t_fin+T_ctrl)/T_ctrl);


% UGV position
figure(2); hold on; grid on;
title('Position');
m = 3; n = 1;

subplot(m, n, 1); hold on; grid on; 
title('x');
plot(t_vec, X(4,1:end), 'k', 'DisplayName', 'x');

subplot(m, n, 2); hold on; grid on; 
title('y');
plot(t_vec, X(5,1:end), 'k', 'DisplayName', 'y');

subplot(m, n, 3); hold on; grid on; 
title('psi');
plot(t_vec, X(6,1:end), 'k', 'DisplayName', 'psi');


% UGV curvilinear
figure(3); 
title('Curvilinear states');
m = 3; n = 1;

subplot(m, n, 1); hold on; grid on;
title('s');
% ylim([x_min(4)-0.5 x_max(4)+0.5]);
plot(t_vec, X_MPC(4,1:end), 'k', 'DisplayName', 's');
% yline(x_min(4), 'r--'); yline(x_max(4), 'r--');

subplot(m, n, 2); hold on; grid on;
title('e_y');
% ylim([x_min(5)+x_min(5)*10/100 x_max(5)+x_max(5)*10/100]);
plot(t_vec, X_MPC(5,1:end), 'k', 'DisplayName', 'e_y');
% yline(x_min(5), 'r--'); yline(x_max(5), 'r--');

subplot(m, n, 3); hold on; grid on;
title('e_{psi}');
% ylim([x_min(6)*180/pi+x_min(6)*180/pi*10/100 x_max(6)*180/pi+x_max(6)*180/pi*10/100]);
plot(t_vec, X_MPC(6,1:end)*180/pi, 'k', 'DisplayName', 'e_{psi}');
% yline(x_min(6)*180/pi, 'r--'); yline(x_max(6)*180/pi, 'r--');


% UGV velocity
figure(4); hold on; grid on;
title('Velocity');
m = 3; n = 1;

subplot(m, n, 1); hold on; grid on;
title('v_x');
% ylim([x_min(1)+x_min(1)*10/100 x_max(1)+x_max(1)*10/100]);
plot(t_vec, X_MPC(1,1:end), 'k', 'DisplayName', 'v_x');
% yline(x_min(1), 'r--'); yline(x_max(1), 'r--');

subplot(m, n, 2); hold on; grid on;
title('v_y');
% ylim([x_min(2)+x_min(2)*10/100 x_max(2)+x_max(2)*10/100]);
plot(t_vec, X_MPC(2,1:end), 'k', 'DisplayName', 'v_y');
% yline(x_min(2), 'r--'); yline(x_max(2), 'r--');

subplot(m, n, 3); hold on; grid on;
title('w_z');
% ylim([x_min(3)*180/pi+x_min(3)*180/pi*10/100 x_max(3)*180/pi+x_max(3)*180/pi*10/100]);
plot(t_vec, X_MPC(3,1:end)*180/pi, 'k', 'DisplayName', 'w_z');
% yline(x_min(3)*180/pi, 'r--'); yline(x_max(3)*180/pi, 'r--');


% UGV input
figure(5); hold on; grid on;
title('Input');
m = 2; n = 1;
subplot(m, n, 1); hold on; grid on;
title('delta_f');
% ylim([u_min(1)*180/pi+u_min(1)*180/pi*10/100 u_max(1)*180/pi+u_max(1)*180/pi*10/100]);
plot(t_vec(1:end-1), U_MPC(1,1:end)*180/pi, 'k', 'DisplayName', 'delta_f');
plot(t_vec(1:end-1), alpha_f(1:end)*180/pi, 'b', 'DisplayName', 'alpha_f');
% yline(u_min(1)*180/pi, 'r--'); yline(u_max(1)*180/pi, 'r--');

subplot(m, n, 2); hold on; grid on;
title('a');
% ylim([u_min(2)+u_min(2)*10/100 u_max(2)+u_max(2)*10/100]);
plot(t_vec(1:end-1), U_MPC(2,1:end), 'k', 'DisplayName', 'a');
% yline(u_min(2), 'r--'); yline(u_max(2), 'r--');

