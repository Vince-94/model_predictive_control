function [mpc_design] = mpc_mtx_fun(mpc_par)

n = mpc_par.n;
p = mpc_par.p;
q = mpc_par.q;
Np = mpc_par.Np;
Nc = mpc_par.Nc;


% prediciton matrices
U0 = zeros(q*Nc, 1);

A_bar = zeros(n*Np, n);
B_bar = zeros(n*Np, q*Nc);
C_bar = zeros(p*Np, n);
D_bar = zeros(p*Np, q*Nc);


for i = 1:Np
    for j = 1:Nc
        A_bar((i-1)*n+1:i*n, 1:n) = mpc_par.A^(i);
        C_bar((i-1)*p+1:i*p, 1:n) = mpc_par.C*mpc_par.A^(i);
        if j <= i
            B_bar((i-1)*n+1:i*n, (j-1)*q+1:j*q) = mpc_par.A^(i-j)*mpc_par.B;
            D_bar((i-1)*p+1:i*p, (j-1)*q+1:j*q) = mpc_par.C*mpc_par.A^(i-j)*mpc_par.B;
        else
            B_bar(i, j) = 0;
            D_bar(i, j) = 0;
        end
    end
end


% objective function matrices
Q_cell = repmat({mpc_par.Q}, Np, 1);
Q_bar = blkdiag(Q_cell{:});
R_cell = repmat({mpc_par.R}, Nc, 1);
R_bar = blkdiag(R_cell{:});


% constraints matrices
X_c = repmat(mpc_par.x_c, Np, 1);
U_c = repmat(mpc_par.u_c, Nc, 1);
dU_c = repmat(mpc_par.du_c, Nc, 1);


% structure
mpc_design.A_bar = A_bar;
mpc_design.B_bar = B_bar;
mpc_design.C_bar = C_bar;
mpc_design.D_bar = D_bar;
mpc_design.Q_bar = Q_bar;
mpc_design.R_bar = R_bar;
mpc_design.X_c = X_c;
mpc_design.U_c = U_c;
mpc_design.dU_c = dU_c;
mpc_design.U0 = U0;
