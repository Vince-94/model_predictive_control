function [mpc] = qp_matrices(prediction_model, Np, Nc, Q, R, S, x_c, u_c, du_c, x0, ref, U_old)

%% Assignment
Ad = prediction_model.Ad;
Bd = prediction_model.Bd;
Kd = prediction_model.Kd;

nx = size(Ad,1);
nu = size(Bd,2);

%% Initialization
A_bar = zeros(nx*Np, nx);
B_bar = zeros(nx*Np, nu*Nc);
K_bar = zeros(nx*Np, 1);
Q_bar = zeros(nx*Np, nx*Np);
R_bar = zeros(nu*Nc, nu*Nc);
X_c = zeros(nx*Np, 2);
U_c = zeros(nu*Nc, 2);
dU_c = zeros(nu*Nc, 2);

Adi = eye(nx);
Kdi = zeros(nx,1);



%% QP
for i = 1:Np
    
    Adi = Ad*Adi;
    Kdi = Kdi + (Ad^(i-1)*Kd);
    
    A_bar((i-1)*nx+1:i*nx, 1:nx) = Adi;
    REF((i-1)*nx+1:i*nx, 1) = ref;
    X_c((i-1)*nx+1:i*nx, 1:2) = x_c;
    Q_bar((i-1)*nx+1:i*nx, (i-1)*nx+1:i*nx) = Q;
    
    for j = 1:Nc
        if j <= i
            B_bar((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = Adi*Bd;
        end
        R_bar((j-1)*nu+1:j*nu, (j-1)*nu+1:j*nu) = R;
        S_bar((j-1)*nu+1:j*nu, (j-1)*nu+1:j*nu) = S;
        U_c((j-1)*nu+1:j*nu, 1:2) = u_c;
        dU_c((j-1)*nu+1:j*nu, 1:2) = du_c;
    end
    
    K_bar((i-1)*nx+1:i*nx, 1) = Kdi;
       
end



%% QP
H = R_bar + B_bar'*Q_bar*B_bar;
mpc.H = (H+H')/2;
mpc.f = (A_bar*x0+K_bar-REF)'*Q_bar*(B_bar);



%% Constraints
% state constraints
A_x = [-B_bar
       +B_bar];
b_x = [-X_c(:,1) + A_bar*x0 + K_bar;
       +X_c(:,2) - A_bar*x0 - K_bar];

% input constraints
A_u = [-eye(nu*Nc)
       +eye(nu*Nc)];
b_u = [-U_c(:,1); 
       +U_c(:,2)];

% input rate constraints
A_du = [-eye(nu*Nc)
        +eye(nu*Nc)];
b_du = [-dU_c(:,1) - U_old; 
        +dU_c(:,2) + U_old];
   

mpc.A_ineq = [A_x; A_u; A_du];
mpc.b_ineq = [b_x; b_u; b_du];

end