function [T,U] = FEM_HeatEq_1D(N, c, f, g, BC_Types, BC_Values,  T_final, num_time_steps, theta)
% Assemble the mass matrix
M = MassMatrix(N+1);

% Determine the length of the time steps
n = N + 1;                       % number of points in space
k = T_final / num_time_steps;    % length of the time steps

% Make a vector of the time partitioning
T = linspace(0, T_final, num_time_steps + 1)';

% Matrix to hold all the solutions
U = zeros(num_time_steps + 1, n);

% Get initial values
x = linspace(0, 1, n);
y = applyIC(g,x);
U(1, :) = y';

% Start time loop
for l = 1:num_time_steps
    % Compute t_theta for the current time step
    t = T(l + 1);
    
    % Get stiffness matrix A and load vector b at time t=t_theta
    A =  StiffnessMatrix(N, @(x) c(x,t));
    b =  LoadVector(N, @(x) f(x,t));
    %apply boundary conditions
    [A, b] = apply_boundary_conditions(A, b, BC_Types, BC_Values,t);

    % Compute the coefficient matrices for the time-stepping scheme
    M_theta = M + theta * k * A;
    M_term = M - (1 - theta) * k * A;
    
    % Compute the right-hand side of Equation (4)
    rhs = M_term * y + k * b;
    
    % Apply boundary conditions
    [M_theta, rhs] = apply_boundary_conditions(M_theta, rhs, BC_Types, BC_Values,t);
    
    % Solve the linear system
    y = M_theta \ rhs;
    
    % Store the solution
    U(l + 1, :) = y';
end
end

