function u_star = adaptive_mpc(x0, ref, Np, Q, R, prediction_model, X_c, U_c)

A = prediction_model.Ad;
B = prediction_model.Bd;
K = prediction_model.Kd;


%% Defining YALMIP decision varaibles
x = sdpvar(size(A,2)*ones(1,Np+1), ones(1,Np+1));
u = sdpvar(size(B,2)*ones(1,Np), ones(1,Np));


%% Cost function
Cost = 0;
for i = 1:Np
    Cost = Cost + (x{i}-ref)'*Q*(x{i}-ref) + u{i}'*R*u{i};
%     Cost = Cost + (x{i}-ref)'*Q*(x{i}-ref) + u{i}'*R*u{i} + (u{i}-u_old)'*S*(u{i}-u_old);
end



%% Constraints
Constraints = [x0 == x{1}];
for i = 1:Np
    Constraints = [Constraints;
                   x{i+1} == A*x{i} + B*u{i} + K;
                   X_c.A * x{i} <= X_c.b + [+K; -K];
                   U_c.A * u{i} <= U_c.b;];
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
u_star = double(u{1});


end