function [SS_local, Q_fun_local, X_local, U_local, J_inf_local, indexSS] = select_local_SS(x, X, U, J_inf, t_lap, numSS_it, numSS_points, D)


%% Select trajectories with the smallest time to reach the goal

% sort the trajectories on the time
[sorted_t_lap, index_t_lap] = sort(t_lap);

if (length(t_lap)<numSS_it)
    numSS_it = length(t_lap);
end

for i = 1:numSS_it
    
    X_sel{i} = X{index_t_lap(i)};
    U_sel{i} = U{index_t_lap(i)};
    J_inf_sel{i} = J_inf{index_t_lap(i)};
    
    one{i} = ones(length(X_sel{i}), 1);
    X_vec{i} = x*one{i}';
    
    diff{i} = D*(X_sel{i} - X_vec{i});
    
    norm = vecnorm(diff{i}, 1, 1);

    [min_norm, index_min_norm] = min(norm);
    
    
    if (index_min_norm - numSS_points/2 > 0)
        indexSSandQfun = [-round(numSS_points/2) + index_min_norm, round(numSS_points/2) + index_min_norm + 1];
    else
        indexSSandQfun = [index_min_norm, index_min_norm + round(numSS_points)];
    end
    
    X_local{i} = X_sel{i}(:,indexSSandQfun(1):indexSSandQfun(2));
    U_local{i} = U_sel{i}(:,indexSSandQfun(1):indexSSandQfun(2));
    J_inf_local{i} = J_inf_sel{i}(:,indexSSandQfun(1):indexSSandQfun(2));
    
    indexSS{i} = indexSSandQfun;
    
end

    
SS_local = [X_local{:}];
Q_fun_local = [J_inf_local{:}];


end