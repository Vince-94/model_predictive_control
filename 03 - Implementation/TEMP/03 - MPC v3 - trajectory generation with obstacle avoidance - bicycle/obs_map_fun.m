function [obs_param] = obs_map_fun(sx, sy)


obs_param.Sx = [sx(1), sy(1), sy(4);
                sx(2), sy(1), sy(3);
                sx(3), sy(1), sy(3);
                sx(4), sy(2), sy(4);
                sx(5), sy(2), sy(4);
                sx(6), sy(1), sy(3);
                sx(7), sy(1), sy(3);
                sx(8), sy(2), sy(4);
                sx(9), sy(2), sy(4);
                sx(10), sy(1), sy(4)];
obs_param.Sy = [sy(1), sx(1), sx(10);
                sy(2), sx(4), sx(5);
                sy(2), sx(8), sx(9);
                sy(3), sx(2), sx(3);
                sy(3), sx(6), sx(7);
                sy(4), sx(1), sx(10)];


