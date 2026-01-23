function [z_c] =  obs_avoidance(z, map)


%% Test values
% clear all
% clc
% z = [1; 1; 0; 0];
% r_act = 2;
% r_safe = 0.25;
% [map] = obs_design(r_act, r_safe);


%% Initialization
x_bounds = [0 20
    0 10];
x = z(1);
y = z(2);
x_c = [NaN NaN];
y_c = [NaN NaN];


%% X_c
k = 1;
for i = 1:length(map.Sx)
    if (abs(x-map.Sx(i,1))>map.r_act)
        continue;
    end
    if (map.Sx(i,2)<=y && y<=map.Sx(i,3))
        x_c(k) = map.Sx(i,1);
        k = k+1;
    end
    if k>2
        break;
    end
end

if (isnan(x_c(1)) && isnan(x_c(2)))
    x_c(:) = [x_bound(1,1) x_bound(1,2)];
elseif (isnan(x_c(2)) && x>x_c(1))
    x_c(1) = x_c(1);
    x_c(2) = x_bound(1,2);
elseif (isnan(x_c(2)) && x<x_c(1))
    x_c(2) = x_c(1);
    x_c(1) = x_bound(1,1);
end




%% Y_c
k = 1;

for i=1:length(map.Sy)
    if (abs(y-map.Sy(i,1))>map.r_act)
        continue;
    end

    if (map.Sy(i,2)<=x && x<=map.Sy(i,3))
        y_c(k) = map.Sy(i,1);
        k = k+1;
    end

    if k>2
        break;
    end
end

if (isnan(y_c(1)) && isnan(y_c(2)))
    y_c(:) = [x_bound(2,1) x_bound(2,2)];
elseif (isnan(y_c(2)) && y>y_c(1))
    y_c(1) = y_c(1);
    y_c(2) = x_bound(2,2);
elseif (isnan(y_c(2)) && y<y_c(1))
    y_c(2) = y_c(1);
    y_c(1) = x_bound(2,1);
end



%% z_c = [x_min x_max; y_min y_max] = [z_min z_max]
z_c = [x_c(1)+map.r_safe, x_c(2)-map.r_safe; 
       y_c(1)+map.r_safe, y_c(2)-map.r_safe];