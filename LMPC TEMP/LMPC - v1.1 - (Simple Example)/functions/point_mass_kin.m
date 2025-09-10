function [A, B, X_c, U_c] = point_mass_kin(Ts)


%% State-space matrices
A = [0 0 1 0; 
    0 0 0 1; 
    0 0 0 0; 
    0 0 0 0];
B = [0 0; 
    0 0;
    1 0; 
    0 1];


% system dimensions
nx = size(A,1);
nu = size(B,2);

% discretization
A = eye(size(A))+A*Ts;
B = B*Ts;




%% System constraitns
x_min = [0
         0
         -5
         -5];
x_max = [20
         10
         5
         5];
u_min = [-10
         -10];
u_max = [10
         10];


% state constraint polihedron
A_x = [-1*eye(nx);
       1*eye(nx)];
b_x = [-x_min; x_max];

X_c = Polyhedron(A_x, b_x);


% input constraint polihedron
A_u = [-1*eye(nu);
       1*eye(nu)];
b_u = [-u_min; u_max];

U_c = Polyhedron(A_u, b_u);


end