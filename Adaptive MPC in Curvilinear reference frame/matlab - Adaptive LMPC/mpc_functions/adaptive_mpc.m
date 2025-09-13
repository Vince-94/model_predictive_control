function [U_star] = adaptive_mpc(x0, ref, U_old, Ts, Np, Nc, Q, R, S, x_c, u_c, du_c, J_inf, SS, PointAndTangent, method)

u_old = U_old(1:2);


%% Prediction
prediction_model = bicycle_pred_model(x0, u_old, Ts, PointAndTangent);


%% Optimization
if method == 0
    mpc = qp_matrices(prediction_model, Np, Nc, Q, R, S, x_c, u_c, du_c, x0, ref, U_old, J_inf, SS);
    options = optimset('Display', 'off');
    U_star = quadprog(mpc.H, mpc.f, mpc.A_ineq, mpc.b_ineq, [], [], [], [], [], options);
elseif method == 1
    U_star = yalmip_mpc(x0, ref, u_old, Np, Q, R, S, prediction_model, x_c, u_c, du_c, J_inf, SS);
    U_star = reshape(U_star, 2*Np, 1);
end


end