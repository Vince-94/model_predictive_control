clear all
close all
clc


%% variables
syms X Y psi 
syms vx vy wz
syms s ey epsi
syms delta_f a
syms lf lr m Iz Cf Cr c


%% equations of motion
% global
X_dot = vx*cos(psi) - vy*sin(psi);
Y_dot = vx*sin(psi) + vy*cos(psi);
psi_dot = wz;

% body
vx_dot = a - 2/m*(Cf*(delta_f - ((vy + lf*wz)/vx))*sin(delta_f)) + wz*vy;
vy_dot = 2/m*(Cf*(delta_f - ((vy + lf*wz)/vx))*cos(delta_f) + Cr*(-((vy - lr*wz)/vx))) - wz*vx;
wz_dot = 2/Iz*(lf*Cf*(delta_f - ((vy + lf*wz)/vx))*cos(delta_f) - lr*Cr*(-((vy - lr*wz)/vx)));

% vx_dot = a - 1/m*(2*Cf*(delta_f - atan2((vy + lf*wz),vx))*sin(delta_f)) + wz*vy;
% vy_dot = 1/m*(2*Cf*(delta_f - atan2((vy + lf*wz),vx))*cos(delta_f) + 2*Cr*(-atan2((vy - lr*wz),vx))) - wz*vx;
% wz_dot = 1/Iz*(lf*2*Cf*(delta_f - atan2((vy + lf*wz),vx))*cos(delta_f) - lr*2*Cr*(-atan2((vy - lr*wz),vx)));


% curvilinear
s_dot = (vx*cos(epsi) - vy*sin(epsi))/(1-c*ey);
ey_dot = vx*sin(epsi) + vy*cos(epsi);
epsi_dot = wz - (vx*cos(epsi) - vy*sin(epsi))/(1-c*ey)*c;



%% global derivatives
dX_X = diff(X_dot, X);
dX_Y = diff(X_dot, Y);
dX_psi = diff(X_dot, psi);
dX_delta_f = diff(X_dot, delta_f);
dX_a = diff(X_dot, a);

dY_X = diff(Y_dot, X);
dY_Y = diff(Y_dot, Y);
dY_psi = diff(Y_dot, psi);
dY_delta_f = diff(Y_dot, delta_f);
dY_a = diff(Y_dot, a);

dpsi_X = diff(psi_dot, X);
dpsi_Y = diff(psi_dot, Y);
dpsi_psi = diff(psi_dot, psi);
dpsi_delta_f = diff(psi_dot, delta_f);
dpsi_a = diff(psi_dot, a);






%% body-global derivatices
dvx_X = diff(vx_dot, X);
dvx_Y = diff(vx_dot, Y);
dvx_psi = diff(vx_dot, psi);

dvy_X = diff(vy_dot, X);
dvy_Y = diff(vy_dot, Y);
dvy_psi = diff(vy_dot, psi);

dwz_X = diff(wz_dot, X);
dwz_Y = diff(wz_dot, Y);
dwz_psi = diff(wz_dot, wz);

dX_vx = diff(X_dot, vx);
dX_vy = diff(X_dot, vy);
dX_wz = diff(X_dot, wz);

dY_vx = diff(Y_dot, vx);
dY_vy = diff(Y_dot, vy);
dY_wz = diff(Y_dot, wz);

dpsi_vx = diff(psi_dot, vx);
dpsi_vy = diff(psi_dot, vy);
dpsi_wz = diff(psi_dot, wz);



%% body-curvilinear derivatices









