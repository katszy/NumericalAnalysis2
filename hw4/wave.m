N = 100; 
T = 1; 
dt = 0.01;

% Define the boundary conditions and initial conditions
BC_Types = ["Neumann", "Neumann"]; 
BC_Values = {@(t) pi*cos(pi*t), @(t) pi*cos(pi - pi*t)}; 
u0 = @(x) sin(pi*x);
ut0 = @(x) -pi*cos(pi*x); 
c = @(x) 1; 
f = @ (x , t ) 0;

h = 1 / N;
x = (0:h:1)';

% Construct the mass matrix and stiffness matrix
M = MassMatrix(N+1);
A = StiffnessMatrix(N, c);

% Apply initial conditions
u = applyIC(u0, x);
ut = applyIC(ut0, x);

%Runge-Kutta 
num_steps = ceil(T / dt);
for k = 1:num_steps
    t = (k - 1) * dt;
  
    [A, b] = apply_boundary_conditions(A, LoadVector(N, f), BC_Types, BC_Values, t);
    
    % Compute the intermediate solution at t + dt/2
    u_mid = u + 0.5 * dt * ut;
    ut_mid = ut + 0.5 * dt * (A * u - M \ b);
    
    % Update the solution at t + dt
    u = u + dt * ut_mid;
    ut = ut + dt * (A * u_mid - M \ b);
end

% Compute the exact solution at the final time
exact_solution = @(x, t) sin(pi*x - pi*t);
u_exact = exact_solution(x, T);

% Plot the exact solution and the approximate solution at the final time
figure;
plot(x, u_exact, 'r--', 'LineWidth', 2);
hold on;
plot(x, u, 'b-', 'LineWidth', 1);
hold off;
legend('Exact Solution', 'Approximate Solution');
xlabel('x');
ylabel('u(x, t)');
title(sprintf('Wave Equation Solution at t = %.2f', T));