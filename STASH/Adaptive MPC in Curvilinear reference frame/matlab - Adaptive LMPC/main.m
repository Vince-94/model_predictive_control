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

%% UGV = [vx, vy, wz, X, Y, psi];
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



%% Track: [x, y, psi, s, length, curvature)
track_n = 4;
W_track = 2;
[PointAndTangent, x0_track] = track_fun(track_n, W_track);




%% ==================================================================
%  ====================== LMPC initialization =======================
%  ==================================================================

%% MPC tuning
method = 1;
T_ctrl = 0.1;
Np = 20;
Nc = Np;
Q = diag([1 1 1 100 1 1]);
R = diag([1 1]);
S = diag([0 0]);
D = diag([0 0 0 1 0 0]);            % select which states take part to build SS_local


%% Feasible solution
load('trajectory.mat')
points = 50;
stack = [X_MPC; X_plant; [[0; 0], U_MPC]];
stack = stack(:, 1:length(stack)/points:end);

X0 = stack(1:6, :);
X0_plant = stack(7:12, :);
U0 = stack(13:14, 2:end);
x0 = X0(:,1);
ref = X0(:,end);
t_lap(1) = t_fin;

figure(1);
plot(X0_plant(4,:), X0_plant(5,:), '-r');


%% Stack initialization
X{1} = X0;
X_global{1} = X0_plant;
U{1} = U0;
J_inf{1} = cost_to_go(X0, U0, ref, Q, R);


% Polyedron
SSQfun = Polyhedron([X{1}', J_inf{1}']);
% SSQfun.computeHRep();
SSQfun.computeVRep();
SS = SSQfun.V(:,1:length(x0))';
Q_fun = SSQfun.V(:, end)';



%% Constraints
% state constraints [vx, vy, wz, s, ey, epsi]
x_c = [0,           +3
       -5,          +5
       -90*pi/180,  +90*pi/180
       0,           +ref(4)+0.5
       -W_track,    W_track
       -90*pi/180,  90*pi/180];

% input constraints [delta_f, a]
u_c = [-25*pi/180,  +25*pi/180
       -1.5,        +1];

% input rate constraints [delta_f, a]
du_c = [-30*pi/180,  +30*pi/180
        -1,          +1];





%% ==================================================================
%  ============================== LMPC ==============================
%  ==================================================================

t_sim = 30;
numSS_it = 4;                               % num of trajectories to be used at each iteration to build the SS
numSS_points = 12;                          % num of near points to select from each trajectory
numSS_points_tot = numSS_points*numSS_it;   % num of points to select from the selected trajectories to build the SS


j = 1;
flag_2 = 0;
while (flag_2 == 0)
    
    
    % iteration initialization
    X_LMPC = x0;
    x_plant = X0_plant(:,1);
    U_old = zeros(Nc*2,1);

    k = 1;
    flag_1 = 0;
    while (flag_1 == 0)
        t = T_ctrl*k;
        clc
        fprintf('Time: t = %.2f, Time step: k = %d, Iteration: j = %d, Initial Iteration Cost: J_iter(0) = %.2f\n', [t, k, j, J_inf{j}(1)]);
        fprintf('Velocity =         [vx = %.2fm/s,  vy = %.2fm/s,   wz =   %.1f°/s] \n', [X_LMPC(1,k), X_LMPC(2,k), X_LMPC(3,k)*180/pi]);
        fprintf('Curvilinear pose = [s =  %.2fm,    ey = %.2fm,     epsi = %.1f°] \n',   [X_LMPC(4,k), X_LMPC(5,k), X_LMPC(6,k)*180/pi]);
        
        
        %% Local convex SS
        [SS_local, Q_fun_local, X_local, U_local, J_inf_local, indexSS] = select_local_SS(X_LMPC(:,k), X, U, J_inf, t_lap, numSS_it, numSS_points, D);
        fprintf('Index SS local = %d \n', [indexSS{1}]);

        %% MPC
        [U_star] = adaptive_mpc(X_LMPC(:,k), ref, U_old, T_ctrl, Np, Nc, Q, R, S, x_c, u_c, du_c, Q_fun_local, SS_local, PointAndTangent, method);
        U_old = U_star;
        if (isempty(U_star) | isnan(U_star))
            u_star = [];
            fprintf('Problem infeasible at k = %d, iteration = %d \n', [k, j]);
            t_fin = T_ctrl*(k-1);
            break;
        end
        U_LMPC(:,k) = U_star(1:2);
        u_old = U_LMPC(:,k);
        fprintf('Input =            [delta_f = %.1f°,  a = %.2fm/s^2]\n', [U_LMPC(1,k)*180/pi, U_LMPC(2,k)]);

        
        
        %% Plant
        x = x_plant(:,k);
        u = U_LMPC(:,k);
        for i=1:T_ctrl/T
            [x_plus, alpha] = bicycle_plant_discrete(x, u, sys_model, T);
            x = x_plus;
            figure(1); grid on; hold on;
            p1(i) = plot(x(4), x(5), '.k');
        end
        x_plant(:,k+1) = x;
        
        
        
        %% Update sample
        k = k+1;

        figure(1); hold on; grid on; axis equal;
        p2 = plot(x_plant(4,k), x_plant(5,k), 'or', 'LineWidth', 0.5);
        arrow(x_plant(4,k), x_plant(5,k), x_plant(6,k), 'b');
        


        %% Conversion from global to curvilinear coordinates
        [s, ey, epsi, CompletedFlag] = global2local_fun(x_plant(4:6,k), x0_track, PointAndTangent, W_track);
        X_LMPC(:,k) = [x_plant(1,k); x_plant(2,k); x_plant(3,k); s; ey; epsi];

        if CompletedFlag == 0
            flag_1 = 1;
            t_fin = T_ctrl*(k);
            break;
        end
        
        
        
        %% Exit condition
        if (abs(X_LMPC(4,k) - ref(4)) <= 10e-2)
            fprintf('Terget reached \n');
            flag_1 = 1;
            t_fin = T_ctrl*(k-1);
        elseif (t == t_sim)
            fprintf('Simulation time reached \n');
            flag_1 = 1;
            t_fin = T_ctrl*(k-1);
        end
        
        

    end
    
    
    
    break;
end















return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


while (flag_2 == 0)


    
    %% Store iteration data
    X{j+1} = X_LMPC;
    X_global{j+1} = X0_plant;
    U{j+1} = U_LMPC;
    J_inf{j+1} = cost_to_go(X_LMPC, U_LMPC, Q, R);
    
    
    %% Update Safe Set and Q-function
    SS = [SS, X_LMPC];
    Q_fun = [Q_fun, J_inf{j+1}];
    
    
    %% Exit condition
    if (abs(J_inf{j}(1) - J_inf{j+1}(1)) < 10e-4)
        flag_2 = 1;
    end
    
    
    
    %% Update iteration
    j = j + 1;
    


    
end
    




%% ===================================================================
%  =========================== Performance ===========================
%  ===================================================================

SS = Polyhedron(SS');

for i = 1:length(X)-1
    
    % Performance
    figure(3); hold on; grid on;
    title('Performance');
    xlabel('$$k$$','interpreter','latex','fontsize',15);
    ylabel('$$J_{\infty}$$','interpreter','latex','fontsize',15);
    plot(find(J_inf{i}(:))-1, J_inf{i}(:), '-');

end



return



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
ylim([x_c(4,1)-0.5, x_c(4,2)+0.5]);
plot(t_vec, X_MPC(4,1:end), 'k', 'DisplayName', 's');
yline(x_c(4,1), 'r--'); yline(x_c(4,2), 'r--');

subplot(m, n, 2); hold on; grid on;
title('e_y');
ylim([x_c(5,1)-x_c(5,1)*10/100, x_c(5,2)+x_c(5,2)*10/100]);
plot(t_vec, X_MPC(5,1:end), 'k', 'DisplayName', 'e_y');
yline(x_c(5,1), 'r--'); yline(x_c(5,2), 'r--');

subplot(m, n, 3); hold on; grid on;
title('e_{psi}');
ylim([x_c(6,1)*180/pi+x_c(6,1)*180/pi*10/100, x_c(6,2)*180/pi+x_c(6,2)*180/pi*10/100]);
plot(t_vec, X_MPC(6,1:end)*180/pi, 'k', 'DisplayName', 'e_{psi}');
yline(x_c(6,1)*180/pi, 'r--'); yline(x_c(6,2)*180/pi, 'r--');


% UGV velocity
figure(4); hold on; grid on;
title('Velocity');
m = 3; n = 1;

subplot(m, n, 1); hold on; grid on;
title('v_x');
ylim([x_c(1,1)+x_c(1,1)*10/100, x_c(1,2)+x_c(1,2)*10/100]);
plot(t_vec, X_MPC(1,1:end), 'k', 'DisplayName', 'v_x');
yline(x_c(1,1), 'r--'); yline(x_c(1,2), 'r--');

subplot(m, n, 2); hold on; grid on;
title('v_y');
ylim([x_c(2,1)+x_c(2,1)*10/100, x_c(2,2)+x_c(2,2)*10/100]);
plot(t_vec, X_MPC(2,1:end), 'k', 'DisplayName', 'v_y');
yline(x_c(2,1), 'r--'); yline(x_c(2,2), 'r--');

subplot(m, n, 3); hold on; grid on;
title('w_z');
ylim([x_c(3,1)*180/pi+x_c(3,1)*180/pi*10/100, x_c(3,2)*180/pi+x_c(3,2)*180/pi*10/100]);
plot(t_vec, X_MPC(3,1:end)*180/pi, 'k', 'DisplayName', 'w_z');
yline(x_c(3,1)*180/pi, 'r--'); yline(x_c(3,2)*180/pi, 'r--');


% UGV input
figure(5); hold on; grid on;
title('Input');
m = 2; n = 1;
subplot(m, n, 1); hold on; grid on;
title('delta_f');
ylim([u_c(1,1)*180/pi+u_c(1,1)*180/pi*10/100, u_c(1,2)*180/pi+u_c(1,2)*180/pi*10/100]);
plot(t_vec(1:end-1), U_MPC(1,1:end)*180/pi, 'k', 'DisplayName', 'delta_f');
plot(t_vec(1:end-1), alpha_f(1:end)*180/pi, 'b', 'DisplayName', 'alpha_f');
yline(u_c(1,1)*180/pi, 'r--'); yline(u_c(1,2)*180/pi, 'r--');

subplot(m, n, 2); hold on; grid on;
title('a');
ylim([u_c(2,1)+u_c(2,1)*10/100 u_c(2,2)+u_c(2,2)*10/100]);
plot(t_vec(1:end-1), U_MPC(2,1:end), 'k', 'DisplayName', 'a');
yline(u_c(2,1), 'r--'); yline(u_c(2,2), 'r--');




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
