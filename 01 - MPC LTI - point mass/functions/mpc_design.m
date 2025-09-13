function [mpc] = mpc_design(mpc_param)

%% Definition
Ts = mpc_param.Ts;
nx = mpc_param.nx;
ny = mpc_param.ny;
nu = mpc_param.nu;
Np = mpc_param.Np;
Nc = mpc_param.Nc;
A = mpc_param.A;
B = mpc_param.B;
C = mpc_param.C;
R = mpc_param.R;
Q = mpc_param.Q;
x_c = mpc_param.x_c;
u_c = mpc_param.u_c;
du_c = mpc_param.du_c;



%% Code

% prediciton matrices
U0 = zeros(nu*Nc, 1);
A_bar = zeros(nx*Np, nx);
B_bar = zeros(nx*Np, nu*Nc);
C_bar = zeros(ny*Np, nx);
D_bar = zeros(ny*Np, nu*Nc);

for i = 1:Np
    for j = 1:Nc
        
        A_bar((i-1)*nx+1:i*nx, 1:nx) = A^(i);
        C_bar((i-1)*ny+1:i*ny, 1:nx) = C*A^(i);
        
        if j <= i
            B_bar((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = A^(i-j)*B;
            D_bar((i-1)*ny+1:i*ny, (j-1)*nu+1:j*nu) = C*A^(i-j)*B;
        else
            B_bar(i, j) = 0;
            D_bar(i, j) = 0;
        end
    end
end


% objective function matrices
Q_cell = repmat({Q}, Np, 1);
Q_bar = blkdiag(Q_cell{:});
R_cell = repmat({R}, Nc, 1);
R_bar = blkdiag(R_cell{:});


% constraints matrices
X_c = repmat(x_c, Np, 1);
U_c = repmat(u_c, Nc, 1);
dU_c = repmat(du_c, Nc, 1);



%% Structure

mpc.Ts = Ts;
mpc.nx = nx;
mpc.ny = ny;
mpc.nu = nu;
mpc.Np = Np;
mpc.Nc = Nc;
mpc.U0 = U0;

mpc.A = A;
mpc.B = B;
mpc.C = C;
mpc.Q = Q;
mpc.R = R;
mpc.A = A;
mpc.x_c = x_c;
mpc.u_c = u_c;
mpc.du_c = du_c;

mpc.A_bar = A_bar;
mpc.B_bar = B_bar;
mpc.C_bar = C_bar;
mpc.D_bar = D_bar;
mpc.Q_bar = Q_bar;
mpc.R_bar = R_bar;
mpc.X_c = X_c;
mpc.U_c = U_c;
mpc.dU_c = dU_c;