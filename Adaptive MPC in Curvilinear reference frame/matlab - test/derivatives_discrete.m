clear all
close all
clc


%% variables
syms X Y psi 
syms vx vy wz
syms s ey epsi
syms delta_f a
syms lf lr m Iz Cf Cr c Ts









%% equations of motion
% vx_plus = vx + (a - 2/m*(Cf*(delta_f - atan2((vy + lf*wz),vx))*sin(delta_f)) + wz*vy)*Ts;
% vy_plus = vy + (2/m*(Cf*(delta_f - atan2((vy + lf*wz),vx))*cos(delta_f) - Cr*atan2((vy - lr*wz),vx)) - wz*vx)*Ts;
% wz_plus = wz + (2/Iz*(lf*Cf*(delta_f - atan2((vy + lf*wz),vx))*cos(delta_f) + lr*Cr*atan2((vy - lr*wz),vx)))*Ts;

vx_plus = vx + (a - 2/m*(Cf*(delta_f - ((vy + lf*wz)/vx))*sin(delta_f)) + wz*vy)*Ts;
vy_plus = vy + (2/m*(Cf*(delta_f - ((vy + lf*wz)/vx))*cos(delta_f) - Cr*((vy - lr*wz)/vx)) - wz*vx)*Ts;
wz_plus = wz + (2/Iz*(lf*Cf*(delta_f - ((vy + lf*wz)/vx))*cos(delta_f) + lr*Cr*((vy - lr*wz)/vx)))*Ts;

s_plus = s + ((vx*cos(epsi) - vy*sin(epsi))/(1-c*ey))*Ts;
ey_plus = ey + (vx*sin(epsi) + vy*cos(epsi))*Ts;
epsi_plus = epsi + (wz - (vx*cos(epsi) - vy*sin(epsi))/(1-c*ey)*c)*Ts;








return

% alpha_f = delta_f - atan2((vy + lf*wz),vx);
% alpha_r = -atan2((vy - lr*wz),vx);
alpha_f = delta_f - ((vy + lf*wz)/vx);
alpha_r = -((vy - lr*wz)/vx);

F_yf = Cf*alpha_f;
F_yr = Cr*alpha_r;


% body
vx_plus = vx + (a - 2/m*(F_yf*sin(delta_f)) + vy*wz)*Ts;
vy_plus = vy + (2/m*(F_yf*cos(delta_f) + F_yr) - wz*vx)*Ts;
wz_plus = wz + (2/Iz*(lf*F_yf*cos(delta_f) - lr*F_yr))*Ts;

% curvilinear
s_plus = s + ((vx*cos(epsi) - vy*sin(epsi))/(1-c*ey))*Ts;
ey_plus = ey + (vx*sin(epsi) + vy*cos(epsi))*Ts;
epsi_plus = epsi + (wz - (vx*cos(epsi) - vy*sin(epsi))/(1-c*ey)*c)*Ts;


%% A11
dvx_vx = diff(vx_plus, vx);
dvx_vy = diff(vx_plus, vy);
dvx_wz = diff(vx_plus, wz);

dvy_vx = diff(vy_plus, vx);
dvy_vy = diff(vy_plus, vy);
dvy_wz = diff(vy_plus, wz);

dwz_vx = diff(wz_plus, vx);
dwz_vy = diff(wz_plus, vy);
dwz_wz = diff(wz_plus, wz);



%% A22
ds_s = diff(s_plus, s);
ds_ey = diff(s_plus, ey);
ds_epsi = diff(s_plus, epsi);

dey_s = diff(ey_plus, s);
dey_ey = diff(ey_plus, ey);
dey_epsi = diff(ey_plus, epsi);

depsi_s = diff(epsi_plus, s);
depsi_ey = diff(epsi_plus, ey);
depsi_epsi = diff(epsi_plus, epsi);




%% A12

dvx_s = diff(vx_plus, s);
dvx_ey = diff(vx_plus, ey);
dvx_epsi = diff(vx_plus, epsi);

dvy_s = diff(vy_plus, s);
dvy_ey = diff(vy_plus, ey);
dvy_epsi = diff(vy_plus, epsi);

dwz_s = diff(wz_plus, s);
dwz_ey = diff(wz_plus, ey);
dwz_epsi = diff(wz_plus, epsi);




%% A21
ds_vx = diff(s_plus, vx);
ds_vy = diff(s_plus, vy);
ds_wz = diff(s_plus, wz);

dey_vx = diff(ey_plus, vx);
dey_vy = diff(ey_plus, vy);
dey_wz = diff(ey_plus, wz);

depsi_vx = diff(epsi_plus, vx);
depsi_vy = diff(epsi_plus, vy);
depsi_wz = diff(epsi_plus, wz);




%% B1
dvx_delta_f = diff(vx_plus, delta_f);
dvx_a = diff(vx_plus, a);

dvy_delta_f = diff(vy_plus, delta_f);
dvy_a = diff(vy_plus, a);

dwz_delta_f = diff(wz_plus, delta_f);
dwz_a = diff(wz_plus, a);

B1 = [dvx_delta_f, dvx_a;
      dvy_delta_f, dvy_a;
      dwz_delta_f, dwz_a];


%% B2
ds_delta_f = diff(s_plus, delta_f);
ds_a = diff(s_plus, a);

dey_delta_f = diff(ey_plus, delta_f);
dey_a = diff(ey_plus, a);

depsi_delta_f = diff(epsi_plus, delta_f);
depsi_a = diff(epsi_plus, a);

B2 = [ds_delta_f,    ds_a;
      dey_delta_f,   dey_a;
      depsi_delta_f, depsi_a];



return




% vx_dot = a - 1/m*(2*Cf*(delta_f - atan2((vy + lf*wz),vx))*sin(delta_f)) + wz*vy;
% vy_dot = 1/m*(2*Cf*(delta_f - atan2((vy + lf*wz),vx))*cos(delta_f) + 2*Cr*(-atan2((vy - lr*wz),vx))) - wz*vx;
% wz_dot = 1/Iz*(lf*2*Cf*(delta_f - atan2((vy + lf*wz),vx))*cos(delta_f) - lr*2*Cr*(-atan2((vy - lr*wz),vx)));




%% body-global derivatices
dvx_X = diff(vx_plus, X);
dvx_Y = diff(vx_plus, Y);
dvx_psi = diff(vx_plus, psi);

dvy_X = diff(vy_plus, X);
dvy_Y = diff(vy_plus, Y);
dvy_psi = diff(vy_plus, psi);

dwz_X = diff(wz_plus, X);
dwz_Y = diff(wz_plus, Y);
dwz_psi = diff(wz_plus, wz);

dX_vx = diff(X_dot, vx);
dX_vy = diff(X_dot, vy);
dX_wz = diff(X_dot, wz);

dY_vx = diff(Y_dot, vx);
dY_vy = diff(Y_dot, vy);
dY_wz = diff(Y_dot, wz);

dpsi_vx = diff(psi_dot, vx);
dpsi_vy = diff(psi_dot, vy);
dpsi_wz = diff(psi_dot, wz);




