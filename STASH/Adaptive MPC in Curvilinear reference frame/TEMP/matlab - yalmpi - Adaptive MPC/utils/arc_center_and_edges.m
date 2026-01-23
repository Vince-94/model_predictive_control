function h = arc_center_and_edges(A, B, C, verso, margin, spec)

% close all
% clear all
% clc
% A = [0; 1];
% B = [3; 4];
% C = [3; 1];
% verso = -1;

R = norm(C-A) + margin;

if verso>0
    a = atan2(A(2)-C(2),A(1)-C(1));
    b = atan2(B(2)-C(2),B(1)-C(1));
    b = mod(b-a,2*pi)+a; % Ensure that arc moves counterclockwise
else
    a = -atan2(A(2)-C(2),A(1)-C(1));
    b = -atan2(B(2)-C(2),B(1)-C(1));
    b = -mod(b-a,2*pi)+a; % Ensure that arc moves counterclockwise
end


t = linspace(a,b,1000);
x = C(1)+R*cos(t);
y = C(2)+R*sin(t);

plot(x, y, spec)
plot(C(1), C(2), 'ko')