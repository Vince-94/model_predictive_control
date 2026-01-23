function pred_model = bicycle_pred_model(x0, u0, Ts, PointAndTangent)

% operative state
vx = x0(1);
vy = x0(2);
wz = x0(3);
s = x0(4);
ey = x0(5);
epsi = x0(6);

% operative input
delta_f = u0(1);
a = u0(2);

% % parameters
% lf = sys.lf;
% lr = sys.lr;
% m = sys.m;
% Iz = sys.Iz;
% Cf = sys.Cf;
% Cr = sys.Cr;


%% Constants
lf = 0.75;
lr = 0.75;
L = lf+lr;
W = 0.8;
m = 800;
Iz = m*(L^2 + W^2)/12;
Cf = 4000;
Cr = 4000;



%% Curvature
index = find(s >= PointAndTangent(:,4) & s < PointAndTangent(:,4)+PointAndTangent(:,5), 1);
c = PointAndTangent(index,6);


%% Operative point

% sideslip angles
alpha_f = delta_f - atan2((vy + lf*wz),vx);
alpha_r = -atan2((vy - lr*wz),vx);


% tire lateral forces
% F_yf = 2*Cf*alpha_f;
% F_yr = 2*Cr*alpha_r;

if (alpha_f > -10*pi/180 && alpha_f < 10*pi/180)
    F_yf = 2*Cf*alpha_f;
else
    alpha_f_lim = 10*pi/180*sign(alpha_f);
    F_yf = 2*Cf*alpha_f_lim;
end

if (alpha_f > -10*pi/180 && alpha_f < 10*pi/180)
    F_yr = 2*Cr*alpha_r;
else
    alphar_r_lim = 10*pi/180*sign(alpha_r);
    F_yr = 2*Cr*alphar_r_lim;
end



% dynamics
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
dvx_vx = -(2*Cf*sin(delta_f)*(vy + lf*wz))/(m*vx^2);
dvx_vy = wz + (2*Cf*sin(delta_f))/(m*vx);
dvx_wz = vy + (2*Cf*lf*sin(delta_f))/(m*vx);

dvy_vx = ((2*Cr*(vy - lr*wz))/vx^2 + (2*Cf*cos(delta_f)*(vy + lf*wz))/vx^2)/m - wz;
dvy_vy = -((2*Cr)/vx + (2*Cf*cos(delta_f))/vx)/m;
dvy_wz = ((2*Cr*lr)/vx - (2*Cf*lf*cos(delta_f))/vx)/m - vx;

dwz_vx = -((2*Cr*lr*(vy - lr*wz))/vx^2 - (2*Cf*lf*cos(delta_f)*(vy + lf*wz))/vx^2)/Iz;
dwz_vy = ((2*Cr*lr)/vx - (2*Cf*lf*cos(delta_f))/vx)/Iz;
dwz_wz = -((2*Cf*cos(delta_f)*lf^2)/vx + (2*Cr*lr^2)/vx)/Iz;                     % WARNING

A11 = [dvx_vx,  dvx_vy,  dvx_wz;
       dvy_vx,  dvy_vy,  dvy_wz;
       dwz_vx,  dwz_vy,  dwz_wz];


   
%% A22
ds_s = 0;
ds_ey = (c*(vx*cos(epsi) - vy*sin(epsi)))/(c*ey - 1)^2;
ds_epsi = (vy*cos(epsi) + vx*sin(epsi))/(c*ey - 1);

dey_s = 0;
dey_ey = 0;
dey_epsi = vx*cos(epsi) - vy*sin(epsi);

depsi_s = 0;
depsi_ey = -(c^2*(vx*cos(epsi) - vy*sin(epsi)))/(c*ey - 1)^2;
depsi_epsi = -(c*(vy*cos(epsi) + vx*sin(epsi)))/(c*ey - 1);

A22 = [ds_s,     ds_ey,     ds_epsi;
       dey_s,    dey_ey,    dey_epsi;
       depsi_s,  depsi_ey,  depsi_epsi];
   
   
%% A12
dvx_s = 0;
dvx_ey = 0;
dvx_epsi = 0;

dvy_s = 0;
dvy_ey = 0;
dvy_epsi = 0;

dwz_s = 0;
dwz_ey = 0;
dwz_epsi = 0;

A12 = [dvx_s,  dvx_ey,  dvx_epsi;
       dvy_s,  dvy_ey,  dvy_epsi;
       dwz_s,  dwz_ey,  dwz_epsi];


%% A21
ds_vx = -cos(epsi)/(c*ey - 1);
ds_vy = sin(epsi)/(c*ey - 1);
ds_wz = 0;

dey_vx = sin(epsi);
dey_vy = cos(epsi);
dey_wz = 0;

depsi_vx = (c*cos(epsi))/(c*ey - 1);
depsi_vy = -(c*sin(epsi))/(c*ey - 1);
depsi_wz = 1;

A21 = [ds_vx,     ds_vy,     ds_wz;
       dey_vx,    dey_vy,    dey_wz;
       depsi_vx,  depsi_vy,  depsi_wz];
   
   
%% B1
dvx_delta_f = - (2*Cf*sin(delta_f))/m - (2*Cf*cos(delta_f)*(delta_f - (vy + lf*wz)/vx))/m;
dvx_a = 1;

dvy_delta_f = (2*Cf*cos(delta_f) - 2*Cf*sin(delta_f)*(delta_f - (vy + lf*wz)/vx))/m;
dvy_a = 0;

dwz_delta_f = (2*Cf*lf*cos(delta_f) - 2*Cf*lf*sin(delta_f)*(delta_f - (vy + lf*wz)/vx))/Iz;
dwz_a = 0;

B1 = [dvx_delta_f,  dvx_a;
      dvy_delta_f,  dvy_a;
      dwz_delta_f,  dwz_a];
  
  
%% B2
ds_delta_f = 0;
ds_a = 0;

dey_delta_f = 0;
dey_a = 0;

depsi_delta_f = 0;
depsi_a = 0;

B2 = [ds_delta_f,    ds_a;
      dey_delta_f,   dey_a;
      depsi_delta_f, depsi_a];

  
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

K = x_dot0 - A*x0 - B*u0;


%% Discretization
Ad = eye(size(A)) + A*Ts;
Bd = B*Ts;
Kd = K*Ts;



%% Structure
pred_model.Ad = Ad;
pred_model.Bd = Bd;
pred_model.Kd = Kd;



end