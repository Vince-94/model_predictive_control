function [mpc] = mpc_design(mpc_param)

%% Definition
Ts = mpc_param.Ts;
n = mpc_param.n;
p = mpc_param.p;
q = mpc_param.q;
Np = mpc_param.Np;
Nc = mpc_param.Nc;
A = mpc_param.A;
B = mpc_param.B;
R = mpc_param.R;
Q = mpc_param.Q;
x_c = mpc_param.x_c;
u_c = mpc_param.u_c;
du_c = mpc_param.du_c;



%% Code

% prediciton matrices
U0 = zeros(q*Nc, 1);
A_bar = zeros(n*Np, n);
B_bar = zeros(n*Np, q*Nc);

for i = 1:Np
    for j = 1:Nc
        
        A_bar((i-1)*n+1:i*n, 1:n) = A^(i);
        %A_bar((i-1)*n+1:i*n, 1:n) = Ad^(i-1);
        
        if j <= i
            B_bar((i-1)*n+1:i*n, (j-1)*q+1:j*q) = A^(i-j)*B;
            %B_bar((i-1)*n+1:i*n, (j-1)*q+1:j*q) = Ad^(i-j-1)*Bd;
        else
            B_bar(i, j) = 0;
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
mpc.n = n;
mpc.p = p;
mpc.q = q;
mpc.Np = Np;
mpc.Nc = Nc;
mpc.U0 = U0;

mpc.A = A;
mpc.B = B;
mpc.Q = Q;
mpc.R = R;
mpc.A = A;
mpc.x_c = x_c;
mpc.u_c = u_c;
mpc.du_c = du_c;

mpc.A_bar = A_bar;
mpc.B_bar = B_bar;
mpc.Q_bar = Q_bar;
mpc.R_bar = R_bar;
mpc.X_c = X_c;
mpc.U_c = U_c;
mpc.dU_c = dU_c;


