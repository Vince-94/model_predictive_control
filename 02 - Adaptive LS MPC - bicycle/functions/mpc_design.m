function [mpc] = mpc_design(mpc_param, pred_model, ref)

%% Assignation
Ts = mpc_param.Ts;
nx = mpc_param.nx;
ny = mpc_param.ny;
nu = mpc_param.nu;
Np = mpc_param.Np;
Nc = mpc_param.Nc;
R = mpc_param.R;
Q = mpc_param.Q;

x_c = mpc_param.x_c;
u_c = mpc_param.u_c;
du_c = mpc_param.du_c;

Ad = pred_model.Ad;
Bd = pred_model.Bd;
Kd = pred_model.Kd;
Cd = pred_model.Cd;
Dd = pred_model.Dd;




%% Initialization
A_bar = zeros(nx*Np, nx);
B_bar = zeros(nx*Np, nu*Nc);
K_bar = zeros(nx*Np, 1);
C_bar = zeros(ny*Np, nx);
D_bar = zeros(ny*Np, nu*Nc);
L_bar = zeros(ny*Np, 1);
Q_bar = zeros(ny*Np, ny*Np);
R_bar = zeros(nu*Nc, nu*Nc);
X_c = zeros(nx*Np, 2);
U_c = zeros(nu*Nc, 2);
dU_c = zeros(nu*Nc, 2);

K1 = Kd;
Ld = Cd*Kd;
L1 = Ld;




%% QP Matrices

for i = 1:Np
    for j = 1:Nc
        
        % State matrices
        A_bar((i-1)*nx+1:i*nx, 1:nx) = Ad^(i);
        C_bar((i-1)*ny+1:i*ny, 1:nx) = Cd*Ad^(i);
        if j <= i
            B_bar((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = Ad^(i-j)*Bd;
            D_bar((i-1)*ny+1:i*ny, (j-1)*nu+1:j*nu) = Cd*Ad^(i-j)*Bd;
        end
        
        % K/L matrice
        if (i>=2)
            Kd = Kd + (Ad^(i-1)*K1);
            Ld = L1 + (Cd*Ad^(i-1)*K1);
        end
        K_bar((i-1)*nx+1:i*nx, 1) = Kd;
        L_bar((i-1)*ny+1:i*ny, 1) = Ld;
        
        % Tuning matrices
        Q_bar((i-1)*ny+1:i*ny, (i-1)*ny+1:i*ny) = Q;
        R_bar((j-1)*nu+1:j*nu, (j-1)*nu+1:j*nu) = R;
        
        % Constraints matrices
        X_c((i-1)*nx+1:i*nx, 1:2) = x_c;
        U_c((j-1)*nu+1:j*nu, 1:2) = u_c;
        dU_c((j-1)*nu+1:j*nu, 1:2) = du_c;
        
    end
end

REF = repmat(ref, Np, 1);


%% Structure
mpc.A_bar = A_bar;
mpc.B_bar = B_bar;
mpc.C_bar = C_bar;
mpc.D_bar = D_bar;
mpc.K_bar = K_bar;
mpc.L_bar = L_bar;

mpc.Q_bar = Q_bar;
mpc.R_bar = R_bar;

mpc.X_min = X_c(:,1);
mpc.X_max = X_c(:,2);
mpc.U_min = U_c(:,1);
mpc.U_max = U_c(:,2);
mpc.dU_min = dU_c(:,1);
mpc.dU_max = dU_c(:,2);

mpc.REF = REF;