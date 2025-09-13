function pred_model = bicycle_pred_model_6states_III(x0, u0, sys, Ts, PointAndTangent)

% operative state
vx = x0(1);
% vx = 1;
vy = x0(2);
wz = x0(3);
s = x0(4);
ey = x0(5);
epsi = x0(6);

% operative input
delta_f = u0(1);
a = u0(2);

% parameters
lf = sys.lf;
lr = sys.lr;
m = sys.m;
Iz = sys.Iz;
Cf = sys.Cf;
Cr = sys.Cr;


%% Curvature
index = find(s >= PointAndTangent(:,4) & s < PointAndTangent(:,4)+PointAndTangent(:,5), 1);
c = PointAndTangent(index,6);


%% Operative point
alpha_f = delta_f - atan2((vy + lf*wz),vx);
alpha_r = -atan2((vy - lr*wz),vx);

F_yf = 2*Cf*alpha_f;
F_yr = 2*Cr*alpha_r;

% continuous
vx_dot = a - 1/m*(F_yf*sin(delta_f)) + wz*vy;
vy_dot = 1/m*(F_yf*cos(delta_f) + F_yr) - wz*vx;
wz_dot = 1/Iz*(lf*F_yf*cos(delta_f) - lr*F_yr);
s_dot = (vx*cos(epsi) - vy*sin(epsi))/(1-c*ey);
ey_dot = vx*sin(epsi) + vy*cos(epsi);
epsi_dot = wz - (vx*cos(epsi) - vy*sin(epsi))/(1-c*ey)*c;

x_dot0 = [vx_dot; vy_dot; wz_dot; s_dot; ey_dot; epsi_dot];


% discrete
vx_plus = vx + Ts*(a - 1/m*F_yf*sin(delta_f) + wz*vy);
vy_plus = vy + Ts*(1/m*(F_yf*cos(delta_f) + F_yr) - wz*vx);
wz_plus = wz + Ts*(1/Iz*(lf*F_yf*cos(delta_f) - lr*F_yr));
s_plus = epsi + Ts*(wz - (vx*cos(epsi) - vy*sin(epsi))/(1-c*ey)*c);
ey_plus = s + Ts*((vx*cos(epsi) - vy*sin(epsi))/(1-c*ey) );
epsi_plus = ey + Ts*(vx*sin(epsi) + vy*cos(epsi));

x_plus0 = [vx_plus; vy_plus; wz_plus; s_plus; ey_plus; epsi_plus];



%% ==================================================================
%  =========================== JACOBIAN =============================
%  ==================================================================


%% A11
dvx_vx = 2*Cf*sin(delta_f)*(vy+lf*wz)/(m*vx^2);
dvx_vy = 2*Cf*sin(delta_f)/(m*vx) + wz;
dvx_wz = 2*Cf*lf*sin(delta_f)/(m*vx) + vy;

dvy_vx = -2/(m*vx^2)*(Cf*(vy+lf*wz)*cos(delta_f) + Cr*(vy-lr*wz)) - wz;
dvy_vy = -2/(m*vx)*(Cf*cos(delta_f) + Cr);
dvy_wz = -2/(m*vx)*(Cf*lf*cos(delta_f) - Cr*lr) - vx;

dwz_vx = 2/(Iz*vx^2)*(Cf*lf*(vy+lf*wz)*cos(delta_f) - Cr*lr*(vy-lr*wz));
dwz_vy = 2/(Iz*vx)*(-Cf*lf*cos(delta_f) + Cr*lr);
dwz_wz = 2/(Iz*vx)*(-Cf*(lf^2)*cos(delta_f) - Cr*lr^2);                     % WARNING

A11 = [dvx_vx,  dvx_vy,  dvx_wz;
       dvy_vx,  dvy_vy,  dvy_wz;
       dwz_vx,  dwz_vy,  dwz_wz];



%% A22
ds_s = 0;
ds_ey = ((vx*cos(epsi) - vy*sin(epsi))/((1-c*ey).^2))*c;
ds_epsi = (-vx*cos(epsi) - vy*sin(epsi))/(1-c*ey);

dey_s = 0;
dey_ey = 0;
dey_epsi = vx*cos(epsi) - vy*sin(epsi);

depsi_s = 0;
depsi_ey = -((c.^2)*(vx*cos(epsi) - vy*sin(epsi)))/((c*ey - 1).^2);
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
ds_vx = cos(epsi)/(1-c*ey);
ds_vy = -sin(epsi)/(1-c*ey);
ds_wz = 0;

dey_vx = sin(epsi);
dey_vy = cos(epsi);
dey_wz = 0;

depsi_vx = -c*cos(epsi)/(1-c*ey);       % WARNING
depsi_vy = c*sin(epsi)/(1-c*ey);
depsi_wz = 1;

A21 = [ds_vx,     ds_vy,     ds_wz;
       dey_vx,    dey_vy,    dey_wz;
       depsi_vx,  depsi_vy,  depsi_wz];



%% B1
dvx_delta_f = -2*Cf/m*(sin(delta_f) + (delta_f - (vy + lf*wz)/vx)*cos(delta_f));
dvx_a = 1;

dvy_delta_f = 2*Cf/m*((-delta_f + (vy+lf*wz)/vx)*sin(delta_f) + cos(delta_f));
% dvy_delta_f = 2*Cf/m*((-delta_f + (vy+lf*wz)/vx)*sin(delta_f));
dvy_a = 0;

dwz_delta_f = 2*Cf*lf/m*((-delta_f + (vy+lf*wz)/vx)*sin(delta_f) + cos(delta_f));
% dwz_delta_f = 2*Cf*lf/m*((-delta_f + (vy+lf*wz)/vx)*sin(delta_f));
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

% Kd = K*Ts;
% Kd = x_plus0 - Ad*x0 - Bd*u0;
Kd = zeros(6,1);


%% Structure
pred_model.Ad = Ad;
pred_model.Bd = Bd;
pred_model.Kd = Kd;


end