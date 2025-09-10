function U_star = yalmip_mpc(x0, ref, u_old, Np, Q, R, S, prediction_model, x_c, u_c, du_c)

A = prediction_model.Ad;
B = prediction_model.Bd;
K = prediction_model.Kd;


%% Defining YALMIP decision varaibles
x = sdpvar(size(A,2)*ones(1,Np+1), ones(1,Np+1));
u = sdpvar(size(B,2)*ones(1,Np), ones(1,Np));


%% Defining YALMIP polytope
A_x = [-1*eye(length(x_c))
       1*eye(length(x_c))];
b_x = [-x_c(:,1);
       x_c(:,2)];
X_c = Polyhedron(A_x, b_x);

A_u = [-1*eye(length(u_c))
       1*eye(length(u_c))];
b_u = [-u_c(:,1); 
       u_c(:,2)];
U_c = Polyhedron(A_u, b_u);

A_du = [-1*eye(length(du_c))
       1*eye(length(du_c))];
b_du = [-du_c(:,1) - u_old; 
       du_c(:,2) + u_old];
dU_c = Polyhedron(A_du, b_du);



%% Cost function
Cost = 0;
for i = 1:Np
%     Cost = Cost + (x{i}-ref)'*Q*(x{i}-ref) + u{i}'*R*u{i};
    Cost = Cost + (x{i}-ref)'*Q*(x{i}-ref) + u{i}'*R*u{i} + (u{i}-u_old)'*S*(u{i}-u_old);
end



%% Constraints
Constraints = [x0 == x{1}];
for i = 1:Np
    Constraints = [Constraints;
                   x{i+1} == A*x{i} + B*u{i} + K;
                   X_c.A * x{i} <= X_c.b;
                   U_c.A * u{i} <= U_c.b
                   dU_c.A * u{i} <= dU_c.b];
end

% for i = 1:Np
%     Constraints = [Constraints;
%                    x{i+1} == A*x{i} + B*u{i} + K;
%                    X_c.A * x{i} <= X_c.b + [+K; -K];
%                    U_c.A * u{i} <= U_c.b;
%                    dU_c.A * u{i} <= dU_c.b + [-u_old; +u_old]];
% end


%% Optimization
options = sdpsettings('verbose',0, 'solver', 'quadprog');
% options.OptimalityTolerance = 1e-15;
% options.StepTolerance = 1e-15;
Problem = optimize(Constraints, Cost, options);
Objective = double(Cost);

U_star = double([u{:}]);
% U_star = double(u{:});


end