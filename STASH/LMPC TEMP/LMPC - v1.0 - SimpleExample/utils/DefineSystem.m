function [A, B, Q, R] = DefineSystem(example)
% Set cost function
Q = diag([1 1 1 1]);            % state cost matrix
R = diag([1 1]);                % input cost matrix


% Define the system matrices
A = [0 0 1 0; 
    0 0 0 1; 
    0 0 0 0; 
    0 0 0 0];
B = [0 0; 
    0 0;
    1 0; 
    0 1];

Ts = 0.1;
A = eye(size(A))+A*Ts;
B = B*Ts;



end