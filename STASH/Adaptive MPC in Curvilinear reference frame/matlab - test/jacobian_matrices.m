% function pred_model = jacobian_matrices(x0, u0, sys, Ts, PointAndTangent)
% % operative state
% vx = x0(1);
% vy = x0(2);
% wz = x0(3);
% s = x0(4);
% ey = x0(5);
% epsi = x0(6);
% 
% % operative input
% delta_f = u0(1);
% a = u0(2);
% 
% % parameters
% lf = sys.lf;
% lr = sys.lr;
% m = sys.m;
% Iz = sys.Iz;
% Cf = sys.Cf;
% Cr = sys.Cr;

close all
clear all
clc

syms vx vy wz s ey epsi delta_f a
syms lf lr m Iz Cf Cr c Ts

x0 = [vx; vy; wz; s; ey; epsi];
u0 = [delta_f; a];


%% Curvature
% index = find(s >= PointAndTangent(:,4) & s < PointAndTangent(:,4)+PointAndTangent(:,5), 1);
% c = PointAndTangent(index,6);


%% Operative point
% alpha_f = delta_f - atan2((vy + lf*wz),vx);
% alpha_r = -atan2((vy - lr*wz),vx);

alpha_f = delta_f - ((vy + lf*wz)/vx);
alpha_r = -((vy - lr*wz)/vx);


F_yf = 2*Cf*alpha_f;
F_yr = 2*Cr*alpha_r;

vx_dot = a - 1/m*(F_yf*sin(delta_f)) + wz*vy;
vy_dot = 1/m*(F_yf*cos(delta_f) + F_yr) - wz*vx;
wz_dot = 1/Iz*(lf*F_yf*cos(delta_f) - lr*F_yr);

s_dot = (vx*cos(epsi) - vy*sin(epsi))/(1-c*ey);
ey_dot = vx*sin(epsi) + vy*cos(epsi);
epsi_dot = wz - (vx*cos(epsi) - vy*sin(epsi))/(1-c*ey)*c;

x_dot0 = [vx_dot; vy_dot; wz_dot; s_dot; ey_dot; epsi_dot];




%% ==================================================================
%  =========================== JACOBIAN =============================
%  ==================================================================

%% A11
dvx_vx = diff(vx_dot, vx);
dvx_vy = diff(vx_dot, vy);
dvx_wz = diff(vx_dot, wz);

dvy_vx = diff(vy_dot, vx);
dvy_vy = diff(vy_dot, vy);
dvy_wz = diff(vy_dot, wz);

dwz_vx = diff(wz_dot, vx);
dwz_vy = diff(wz_dot, vy);
dwz_wz = diff(wz_dot, wz);

A11 = [dvx_vx,  dvx_vy,  dvx_wz;
       dvy_vx,  dvy_vy,  dvy_wz;
       dwz_vx,  dwz_vy,  dwz_wz];


%% A22
ds_s = diff(s_dot, s);
ds_ey = diff(s_dot, ey);
ds_epsi = diff(s_dot, epsi);

dey_s = diff(ey_dot, s);
dey_ey = diff(ey_dot, ey);
dey_epsi = diff(ey_dot, epsi);

depsi_s = diff(epsi_dot, s);
depsi_ey = diff(epsi_dot, ey);
depsi_epsi = diff(epsi_dot, epsi);

A22 = [ds_s,     ds_ey,     ds_epsi;
       dey_s,    dey_ey,    dey_epsi;
       depsi_s,  depsi_ey,  depsi_epsi];



%% A12
dvx_s = diff(vx_dot, s);
dvx_ey = diff(vx_dot, ey);
dvx_epsi = diff(vx_dot, epsi);

dvy_s = diff(vy_dot, s);
dvy_ey = diff(vy_dot, ey);
dvy_epsi = diff(vy_dot, epsi);

dwz_s = diff(wz_dot, s);
dwz_ey = diff(wz_dot, ey);
dwz_epsi = diff(wz_dot, epsi);

A12 = [dvx_s,  dvx_ey,  dvx_epsi;
       dvy_s,  dvy_ey,  dvy_epsi;
       dwz_s,  dwz_ey,  dwz_epsi];
   
   
%% A21
ds_vx = diff(s_dot, vx);
ds_vy = diff(s_dot, vy);
ds_wz = diff(s_dot, wz);

dey_vx = diff(ey_dot, vx);
dey_vy = diff(ey_dot, vy);
dey_wz = diff(ey_dot, wz);

depsi_vx = diff(epsi_dot, vx);
depsi_vy = diff(epsi_dot, vy);
depsi_wz = diff(epsi_dot, wz);

A21 = [ds_vx,     ds_vy,     ds_wz;
       dey_vx,    dey_vy,    dey_wz;
       depsi_vx,  depsi_vy,  depsi_wz];
   

%% B1
dvx_a = diff(vx_dot, a);
dvx_delta_f = diff(vx_dot, delta_f);

dvy_a = diff(vy_dot, a);
dvy_delta_f = diff(vy_dot, delta_f);

dwz_a = diff(wz_dot, a);
dwz_delta_f = diff(wz_dot, delta_f);

B1 = [dvx_delta_f,  dvx_a;
      dvy_delta_f,  dvy_a;
      dwz_delta_f,  dwz_a];


%% B2
ds_a = diff(s_dot, a);
ds_delta_f = diff(s_dot, delta_f);

dey_a = diff(ey_dot, a);
dey_delta_f = diff(ey_dot, delta_f);

depsi_a = diff(epsi_dot, a);
depsi_delta_f = diff(epsi_dot, delta_f);

B2 = [ds_delta_f,    ds_a;
      dey_delta_f,   dey_a;
      depsi_delta_f, depsi_a];




%% ==================================================================
%  =========================== QP =============================
%  ==================================================================

%% QP matrices
A = [A11 A12
      A21 A22];
A(isnan(A)) = 0;
A(A==+inf) = 0;
A(A==-inf) = 0;

B = [B1;
      B2];
B(isnan(B)) = 0;
B(B==+inf) = 0;
B(B==-inf) = 0;

K = x_dot0 + A*x0 + B*u0;



%% Discretization
Ad = eye(size(A)) + A*Ts;
Bd = B*Ts;
Kd = K*Ts;



%% Structure
% pred_model.A = A;
% pred_model.B = B;
% pred_model.K = K;
pred_model.Ad = Ad;
pred_model.Bd = Bd;
pred_model.Kd = Kd;


%%
z = sin(delta_f);
zz = diff(z, delta_f)


% end