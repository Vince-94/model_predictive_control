function [mpc_design] = mpc_mtx_fun(mpc_par)


n = mpc_par.n;
p = mpc_par.p;
q = mpc_par.q;
Np = mpc_par.Np;


% prediciton matrices
U0 = zeros(q*Np, 1);

A_bar = zeros(n*(Np+1), n);
B_bar = zeros(n*(Np+1), q*Np);

for i = 1:Np+1
    for j = 1:Np
        
        A_bar((i-1)*n+1:i*n, 1:n) = mpc_par.A^(i-1);
        
        if j < i
            B_bar((i-1)*n+1:i*n, (j-1)*q+1:j*q) = mpc_par.A^(i-j-1)*mpc_par.B;
        else
            B_bar(i,j) = 0;
        end
    end
end


% objective function matrices
Q_cell = repmat({mpc_par.Q}, Np+1, 1);
Q_bar = blkdiag(Q_cell{:});
R_cell = repmat({mpc_par.R}, Np, 1);
R_bar = blkdiag(R_cell{:});


% constraints matrices
X_c = repmat(mpc_par.x_c, Np+1, 1);
U_c = repmat(mpc_par.u_c, Np, 1);
dU_c = repmat(mpc_par.du_c, Np, 1);


% structure
mpc_design.A_bar = A_bar;
mpc_design.B_bar = B_bar;
mpc_design.Q_bar = Q_bar;
mpc_design.R_bar = R_bar;
mpc_design.X_c = X_c;
mpc_design.U_c = U_c;
mpc_design.dU_c = dU_c;
mpc_design.U0 = U0;




