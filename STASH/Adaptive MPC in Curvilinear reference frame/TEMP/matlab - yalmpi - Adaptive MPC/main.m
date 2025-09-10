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
x0 = [2; 0; 0; 0; 0; 0];
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

track_n = 4;
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
x_min = [0; -0.3; -26*pi/180; 0; -W_track; -10*pi/180];
x_max = [3; +0.3; +26*pi/180; +s_target+0.5; +W_track; +10*pi/180];


% input constraints [delta_f, a]
u_min = [-20*pi/180; -0.5];
u_max = [+20*pi/180; +0.2];


% input rate constraints [ddelta_f, da]
du_min = [-10e5; -10e5];
du_max = [+10e5; +10e5];

% Prediction model
prediction_model = bicycle_pred_model(X_MPC, u0, T_ctrl, PointAndTangent);

% constraints polytope
A_x = [-1*eye(length(x_min))
       1*eye(length(x_max))];
b_x = [-x_min;
       x_max];
X_c = Polyhedron(A_x, b_x);

A_u = [-1*eye(length(u_min))
       1*eye(length(u_max))];
b_u = [-u_min; 
       u_max];
U_c = Polyhedron(A_u, b_u);

A_du = [-1*eye(length(du_min))
       1*eye(length(du_max))];
b_du = [-du_min; 
       du_max];
dU_c = Polyhedron(A_du, b_du);




%% ==================================================================
%  ============================== MPC ===============================
%  ==================================================================

% MPC tuning
Np = 15;
Nc = Np;
Q = diag([1 10 10 1 10 10]);
R = diag([100 1]);
S = diag([0 0]);


k = 1;
flag = 0;
t_sim = 15;


while (flag == 0)
    
    clc
    t = T_ctrl*k;

    %% Print info
    fprintf('Sample: k = %d, time: t = %.2f\n', [k, t]);
    fprintf('Np = %d, Nc = %d \n', [Np, Nc]);
    fprintf('Velocity =         [vx = %.2fm/s,  vy = %.2fm/s,   wz =   %.1f°/s] \n', [X_MPC(1,k), X_MPC(2,k), X_MPC(3,k)*180/pi]);
    fprintf('Curvilinear pose = [s =  %.2fm,    ey = %.2fm,     epsi = %.1f°] \n',   [X_MPC(4,k), X_MPC(5,k), X_MPC(6,k)*180/pi]);
%     fprintf('Global pose =      [x =  %.2fm,    y =  %.2fm,     psi =  %.1f°] \n',   [X(4,k), X(5,k), X(6,k)*180/pi]);

    
    %% Optimization
%     U_MPC(:,k) = [20*pi/180; 0.2];
    U_MPC(:,k) = adaptive_mpc(X_MPC(:,k), ref, Np, Q, R, prediction_model, X_c, U_c);
    fprintf('Input =            [delta_f = %.1f°,  a = %.2fm/s^2]\n', [U_MPC(1,k)*180/pi, U_MPC(2,k)]);
    
    if isnan(U_MPC(:,k))
        fprintf('Problem infeasible at k = %d \n', k);
        U_MPC(:,k) = [];
%         t_fin = t + T_ctrl;
        t_fin = T_ctrl*(k-1);
        break;
    end



    %% Plant discrete
    x = X(:,k);
    u = U_MPC(:,k);
    for i=1:T_ctrl/T
        [x_plus, alpha] = bicycle_plant_discrete(x, u, sys_model, T);
        x = x_plus;
%         figure(1); grid on; hold on;
%         p1(i) = plot(x(4), x(5), '.k');
    end
    
    X(:,k+1) = x;
    alpha_f(k) = alpha;

    
    %% Update sample
    k = k+1;
    u_old = U_MPC(:,k-1);
    
    figure(1); hold on; grid on; axis equal;
    p2 = plot(X(4,k), X(5,k), 'or', 'LineWidth', 0.5);
    arrow(X(4,k), X(5,k), X(6,k), 'b');
    
%     delete(p1)
    
    
    %% Conversion from global to curvilinear coordinates
    [s, ey, epsi, CompletedFlag] = global2local_fun(X(4:6,k), x0_track, PointAndTangent, W_track);
    X_MPC(:,k) = [X(1,k); X(2,k); X(3,k); s; ey; epsi];
    
    if CompletedFlag == 0
        flag = 1;
        t_fin = T_ctrl*(k);
        break;
    end
    
    
    %% Prediction
    prediction_model = bicycle_pred_model(X_MPC(:,k), U_MPC(:,k-1), T_ctrl, PointAndTangent);
    
    
    %% Exit condition
    if (abs(X_MPC(4,k) - s_target) <= 10e-2)
        fprintf('Terget reached \n');
        flag = 1;
        t_fin = T_ctrl*(k-1);
    elseif (t == t_sim)
        fprintf('Simulation time reached \n');
        flag = 1;
        t_fin = T_ctrl*(k-1);
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
ylim([x_min(4)-0.5 x_max(4)+0.5]);
plot(t_vec, X_MPC(4,1:end), 'k', 'DisplayName', 's');
yline(x_min(4), 'r--'); yline(x_max(4), 'r--');

subplot(m, n, 2); hold on; grid on;
title('e_y');
ylim([x_min(5)+x_min(5)*10/100 x_max(5)+x_max(5)*10/100]);
plot(t_vec, X_MPC(5,1:end), 'k', 'DisplayName', 'e_y');
yline(x_min(5), 'r--'); yline(x_max(5), 'r--');

subplot(m, n, 3); hold on; grid on;
title('e_{psi}');
ylim([x_min(6)*180/pi+x_min(6)*180/pi*10/100 x_max(6)*180/pi+x_max(6)*180/pi*10/100]);
plot(t_vec, X_MPC(6,1:end)*180/pi, 'k', 'DisplayName', 'e_{psi}');
yline(x_min(6)*180/pi, 'r--'); yline(x_max(6)*180/pi, 'r--');


% UGV velocity
figure(4); hold on; grid on;
title('Velocity');
m = 3; n = 1;

subplot(m, n, 1); hold on; grid on;
title('v_x');
ylim([x_min(1)+x_min(1)*10/100 x_max(1)+x_max(1)*10/100]);
plot(t_vec, X_MPC(1,1:end), 'k', 'DisplayName', 'v_x');
yline(x_min(1), 'r--'); yline(x_max(1), 'r--');

subplot(m, n, 2); hold on; grid on;
title('v_y');
ylim([x_min(2)+x_min(2)*10/100 x_max(2)+x_max(2)*10/100]);
plot(t_vec, X_MPC(2,1:end), 'k', 'DisplayName', 'v_y');
yline(x_min(2), 'r--'); yline(x_max(2), 'r--');

subplot(m, n, 3); hold on; grid on;
title('w_z');
ylim([x_min(3)*180/pi+x_min(3)*180/pi*10/100 x_max(3)*180/pi+x_max(3)*180/pi*10/100]);
plot(t_vec, X_MPC(3,1:end)*180/pi, 'k', 'DisplayName', 'w_z');
yline(x_min(3)*180/pi, 'r--'); yline(x_max(3)*180/pi, 'r--');


% UGV input
figure(5); hold on; grid on;
title('Input');
m = 2; n = 1;
subplot(m, n, 1); hold on; grid on;
title('delta_f');
ylim([u_min(1)*180/pi+u_min(1)*180/pi*10/100 u_max(1)*180/pi+u_max(1)*180/pi*10/100]);
plot(t_vec(1:end-1), U_MPC(1,1:end)*180/pi, 'k', 'DisplayName', 'delta_f');
plot(t_vec(1:end-1), alpha_f(1:end)*180/pi, 'b', 'DisplayName', 'alpha_f');
yline(u_min(1)*180/pi, 'r--'); yline(u_max(1)*180/pi, 'r--');

subplot(m, n, 2); hold on; grid on;
title('a');
ylim([u_min(2)+u_min(2)*10/100 u_max(2)+u_max(2)*10/100]);
plot(t_vec(1:end-1), U_MPC(2,1:end), 'k', 'DisplayName', 'a');
yline(u_min(2), 'r--'); yline(u_max(2), 'r--');




%% Save plots
prompt = 'Do you want to save the plots? y/n [n]: ';
str = input(prompt,'s');
if str == 'y'
    % mpc params
    mpc.Np = Np;
    mpc.Nc = Nc;
    mpc.Q = Q;
    mpc.R = R;
    save data/mpc_data.mat mpc sys_model
    
    % trajectory plot
    filename = sprintf('data/trajectory.png');
    exportgraphics(figure(1), filename);
    
    % trajectory plot
    filename = sprintf('data/velocities.png');
    exportgraphics(figure(3), filename);
    
    % trajectory plot
    filename = sprintf('data/curvliniear_coordinates.png');
    exportgraphics(figure(4), filename);
    
    % trajectory plot
    filename = sprintf('data/input.png');
    exportgraphics(figure(5), filename);
    
end

close all
