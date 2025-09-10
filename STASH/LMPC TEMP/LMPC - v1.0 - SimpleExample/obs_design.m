function [map] = obs_design(r_act, r_safe)

% obstacles
sx = [0; 2; 4; 6; 8; 10; 12; 14; 16; 20];
sy = [0; 3; 7; 10];

%% Map parametrization
map.Sx = [sx(1), sy(1), sy(4);
        sx(2), sy(1), sy(3);
        sx(3), sy(1), sy(3);
        sx(4), sy(2), sy(4);
        sx(5), sy(2), sy(4);
        sx(6), sy(1), sy(3);
        sx(7), sy(1), sy(3);
        sx(8), sy(2), sy(4);
        sx(9), sy(2), sy(4);
        sx(10), sy(1), sy(4)];
map.Sy = [sy(1), sx(1), sx(10);
        sy(2), sx(4), sx(5);
        sy(2), sx(8), sx(9);
        sy(3), sx(2), sx(3);
        sy(3), sx(6), sx(7);
        sy(4), sx(1), sx(10)];
    
map.r_act = r_act;
map.r_safe = r_safe;



