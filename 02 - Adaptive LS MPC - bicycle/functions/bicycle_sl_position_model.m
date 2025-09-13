function [pred_model] = bicycle_sl_position_model(Ts, v_des, p)

Lr = 0.75;

theta_0 = p(1);
delta_0 = p(2);


% Model
Ad = [1 0 0
      0 1 0
      0 0 1];
Bd = [-v_des*sin(theta_0+delta_0)*Ts
      v_des*cos(theta_0+delta_0)*Ts
      v_des/Lr*cos(delta_0)*Ts];

Kd = [v_des*cos(theta_0+delta_0)*Ts
      v_des*sin(theta_0+delta_0)*Ts
      v_des*cos(delta_0)*Ts];

Cd = [1 0 0
    0 1 0];
Dd = zeros(size(Cd,1), size(Bd,2));




%% Structure

pred_model.Ad = Ad;
pred_model.Bd = Bd;
pred_model.Kd = Kd;
pred_model.Cd = Cd;
pred_model.Dd = Dd;