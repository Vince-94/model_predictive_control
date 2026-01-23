function [x_dot] = bicycle_plant_continuous(x, u, sys)

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
lf = sys.lf;
lr = sys.lr;
m = sys.m;
Iz = sys.Iz;
Cf = sys.Cf;
Cr = sys.Cr;



%% Linear tire's model
alpha_f = delta_f - atan2((vy + lf*wz),vx);
alpha_r = -atan2((vy - lr*wz),vx);
% alpha_f = atan2((vy + lf*wz - delta_f*vx),vx);
% alpha_r = atan2((vy - lr*wz),vx);

alpha_lim = 10*pi/180;

if (alpha_f <= alpha_lim && alpha_f >= -alpha_lim)
    F_yf = 2*Cf*alpha_f;
else
    F_yf = 2*Cf*alpha_lim*sign(alpha_f);
end

if (alpha_r <= alpha_lim && alpha_r >= -alpha_lim)
    F_yr = 2*Cr*alpha_r;
else
    F_yr = 2*Cr*alpha_lim*sign(alpha_r);
end



%% State equation
vx_dot = a - 1/m*(F_yf*sin(delta_f)) + wz*vy;
vy_dot = 1/m*(F_yf*cos(delta_f) + F_yr) - wz*vx;
wz_dot = 1/Iz*(lf*F_yf*cos(delta_f) - lr*F_yr);

X_dot = vx*cos(psi) - vy*sin(psi);
Y_dot = vx*sin(psi) + vy*cos(psi);
psi_dot = wz;

x_dot = [vx_dot; vy_dot; wz_dot; X_dot; Y_dot; psi_dot];