function [x_dot, x_plus] = bicycle_plant(x, u, Ts, sys_model)

% state
vx = x(1);
vy = x(2);
wz = x(3);
X = x(4);
Y = x(5);
psi = x(6);

% input
delta_f = u(1);
a = u(2);

% system constants
lf = sys_model.lf;
lr = sys_model.lr;
m = sys_model.m;
Iz = sys_model.Iz;
Cf = sys_model.Cf;
Cr = sys_model.Cr;



%% Linear tire's model
alpha_f = delta_f - atan2((vy + lf*wz),vx);
alpha_r = -atan2((vy - lr*wz),vx);

F_yf = 2*Cf*alpha_f;
F_yr = 2*Cr*alpha_r;


%% Continue plant
% state equation
vx_dot = a - 1/m*(F_yf*sin(delta_f)) + wz*vy;
vy_dot = 1/m*(F_yf*cos(delta_f) + F_yr) - wz*vx;
wz_dot = 1/Iz*(lf*F_yf*cos(delta_f) - lr*F_yr);

X_dot = vx*cos(psi) - vy*sin(psi);
Y_dot = vx*sin(psi) + vy*cos(psi);
psi_dot = wz;

x_dot = [vx_dot; vy_dot; wz_dot; X_dot; Y_dot; psi_dot];



%% Discrete plant
vx_plus = vx + (a - 1/m*(F_yf*sin(delta_f)) + wz*vy)*Ts;
vy_plus = vy + (1/m*(F_yf*cos(delta_f) + F_yr) - wz*vx)*Ts;
wz_plus = wz + (1/Iz*(lf*F_yf*cos(delta_f) - lr*F_yr))*Ts;

X_plus = X + (vx*cos(psi) - vy*sin(psi))*Ts;
Y_plus = Y + (vx*sin(psi) + vy*cos(psi))*Ts;
psi_plus = psi + wz*Ts;

x_plus = [vx_plus; vy_plus; wz_plus; X_plus; Y_plus; psi_plus];