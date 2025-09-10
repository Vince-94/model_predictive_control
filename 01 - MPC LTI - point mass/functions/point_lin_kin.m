function [pred_model] = point_lin_kin(Ts)

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];
C = eye(4);
D = zeros(size(C,1), size(B,2));

% Discretization
Ad = eye(size(A))+A*Ts;
Bd = B*Ts;
Cd = C;
Dd = D;


%% Costraints
x_bound = [-inf inf         % x bounds [m]
           -inf inf         % y bounds [m]
           -1 1             % vx bounds [m/s]
           -1 1];           % vy bounds [m/s]
u_bound = [-0.5 0.5         % ax bounds [m/s^2]
           -0.5 0.5];       % ay bounds [m/s^2]
du_bound = [-1 1
            -1 1];



%% Structure
pred_model.A = A;
pred_model.B = B;
pred_model.C = C;
pred_model.D = D;

pred_model.Ad = Ad;
pred_model.Bd = Bd;
pred_model.Cd = Cd;
pred_model.Dd = Dd;

pred_model.x_bound = x_bound;
pred_model.u_bound = u_bound;
pred_model.du_bound = du_bound;