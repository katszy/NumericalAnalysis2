function [T, U] = SolveTheta(N, c, f, g, BC_Types, BC_Values,  T_final, num_time_steps, theta)
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
    % This is done by interpolating the initial value function g
    x = linspace(0, 1, n);
    u = zeros(n, 1);
    for i = 1:n
        u(i) = g(x(i));
    end
    U(1, :) = u';

    % Start time loop
    for l = 1:num_time_steps
        % Compute t_theta for the current time step
        t = T(l + 1);

        % Get stiffness matrix A and load vector b at time t=t_theta
        A =  StiffnessMatrix(N, @(x) c(x,t));
        b =  LoadVector(N, @(x) f(x,t));
        
        % Account for boundary conditions
            if BC_Types(1) == 'N'
                b(1) = b(1) - BC_Values{1}(t);
            end
            if BC_Types(2) == 'N'
                b(end) = b(end) + BC_Values{2}(t);
            end 

        % Compute the coefficient matrices for the time-stepping scheme
        M_theta = M + theta * k * A;
        M_term = M - (1 - theta) * k * A;

        % Compute the right-hand side of Equation (4)
        rhs = M_term * u + k * b;

        % Apply boundary conditions
        if BC_Types(1) == 'D'
            M_theta(1, :) = 0;
            M_theta(1, 1) = 1;
            rhs(1) = BC_Values{1}(t);
        end
        if BC_Types(2) == 'D'
            M_theta(n, :) = 0;
            M_theta(n, n) = 1;
            rhs(n) = BC_Values{2}(t);
        end

        % Solve the linear system
        u = M_theta \ rhs;

        % Store the solution
        U(l + 1, :) = u';
    end
end
