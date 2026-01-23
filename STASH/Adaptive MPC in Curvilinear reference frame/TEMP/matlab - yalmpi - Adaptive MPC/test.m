%% MAIN
close all
clear all
clc


Np = 20;
Nc = 10;
mpc.m = 3;
mpc.lf = 0.75;


% mpc params
save example.mat Np Nc mpc

return
LastName = {'Sanchez';'Johnson';'Li';'Diaz';'Brown'};
Age = [38;43;38;40;49];
Smoker = logical([1;0;1;0;1]);
Height = [71;69;64;67;64];
Weight = [176;163;131;133;119];
BloodPressure = [124 93; 109 77; 125 83; 117 75; 122 80];

T = table(LastName,Age,Smoker,Height,Weight,BloodPressure);



% fig = uifigure('Name','MPC Parameters');
% fig.WindowStyle = 'modal';
% uit = uitable(fig,'Data',T);
% uit.Position = [20 20 520 320];
% uit.RowName = 'numbered';



return

%%
close all
clear all
clc

myData = magic(100);

imwrite(uint8(myData), 'my image.png')
