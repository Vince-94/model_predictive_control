function [x_plus, alpha_f] = bicycle_plant_discrete(x, u, sys, Ts)

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

theta_f = atan2((vy + lf*wz),vx);
theta_r = atan2((vy - lr*wz),vx);

alpha_f = delta_f - theta_f;
alpha_r = -theta_r;


% alpha_lim = 10*pi/180;
% if (alpha_f <= alpha_lim && alpha_f >= -alpha_lim)
%     F_yf = 2*Cf*alpha_f;
% else
%     F_yf = 2*Cf*alpha_lim*sign(alpha_f);
% end
% 
% if (alpha_r <= alpha_lim && alpha_r >= -alpha_lim)
%     F_yr = 2*Cr*alpha_r;
% else
%     F_yr = 2*Cr*alpha_lim*sign(alpha_r);
% end

F_yf = 2*Cf*alpha_f;
F_yr = 2*Cr*alpha_r;


%% State equation
vx_plus = vx + (a - 1/m*(F_yf*sin(delta_f)) + wz*vy)*Ts;
vy_plus = vy + (1/m*(F_yf*cos(delta_f) + F_yr) - wz*vx)*Ts;
wz_plus = wz + (1/Iz*(lf*F_yf*cos(delta_f) - lr*F_yr))*Ts;

X_plus = X + (vx*cos(psi) - vy*sin(psi))*Ts;
Y_plus = Y + (vx*sin(psi) + vy*cos(psi))*Ts;
psi_plus = psi + wz*Ts;

x_plus = [vx_plus; vy_plus; wz_plus; X_plus; Y_plus; psi_plus];