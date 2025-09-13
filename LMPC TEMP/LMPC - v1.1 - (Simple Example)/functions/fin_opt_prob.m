function [u_star] = fin_opt_prob(x0, xf, Np, Q, R, J_inf, SS, A, B, X_c, U_c, map)


%% Defining YALMIP decision varaibles
x = sdpvar(size(A,2)*ones(1,Np+1), ones(1,Np+1));
u = sdpvar(size(B,2)*ones(1,Np), ones(1,Np));
lambda = sdpvar(length(J_inf), 1);



%% Cost function

% Stage Cost
J_stage = 0;
for i = 1:Np
    J_stage = J_stage + x{i}'*Q*x{i} + u{i}'*R*u{i};
%     J_stage = J_stage + (x{i}-xf)'*Q*(x{i}-xf) + u{i}'*R*u{i};
end

% Complexive cost
Cost = J_stage + J_inf*lambda;



%% Constraints

% System Dynamics Constraints
Constraints = [x0 == x{1}];

for i = 1:Np
    
%     [z_c] = obs_avoidance(x0, x_bound, map);
%     z_min = z_c(:,1);
%     z_max = z_c(:,2);
%     X_c.b(1:2) = z_min;
%     X_c.b(5:6) = z_max;
    
    Constraints = [Constraints;
                 x{i+1} == A*x{i} + B*u{i};
                 X_c.A * x{i} <= X_c.b;
                 U_c.A * u{i} <= U_c.b];
end

% Terminal Constraint: enforce predicted state in the convex safe set
Constraints=[Constraints;
             lambda >= 0;                        % Constraint the multipliers to be positive
             x{Np+1} == SS*lambda;                % Terminal point in the convex hull
             ones(1,length(J_inf))*lambda == 1];  % Must be convex combination --> sum to 1


%% Optimization
options = sdpsettings('verbose',0, 'solver', 'quadprog');
% options.OptimalityTolerance = 1e-15;
% options.StepTolerance = 1e-15;
Problem = optimize(Constraints, Cost, options);
Objective = double(Cost);
u_star = double(u{1});

end