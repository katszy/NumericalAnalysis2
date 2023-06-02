function [T,U] = FEM_HeatEq_1D(N, c, f, g, BC_Types, BC_Values,  T_final, num_time_steps, theta)
M = MassMatrix(N+1);
n = N + 1;                       % number of points in space
k = T_final / num_time_steps;    % length of the time steps

% Time partitioning
T = linspace(0, T_final, num_time_steps + 1)';
% Matrix to hold all the solutions
U = zeros(num_time_steps + 1, n);

% Initial values
x = linspace(0, 1, n);
y = applyIC(g,x);
U(1, :) = y';

% Start time loop
for l = 1:num_time_steps
    t = T(l + 1);
    
    A =  StiffnessMatrix(N, @(x) c(x,t));
    b =  LoadVector(N, @(x) f(x,t));
    [A, b] = apply_boundary_conditions(A, b, BC_Types, BC_Values,t);
    
    M_theta = M + theta * k * A;
    M_term = M - (1 - theta) * k * A;
    rhs = M_term * y + k * b;
   
    [M_theta, rhs] = apply_boundary_conditions(M_theta, rhs, BC_Types, BC_Values,t);
   
    y = M_theta \ rhs;
    U(l + 1, :) = y';
end
end

