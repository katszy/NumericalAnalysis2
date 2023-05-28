% Parameters
N = 100;                  % Number of elements
T = 1;                    % Final time
dt = 0.01;                % Time step size
c = @(x) 1;

% Define the boundary conditions and initial conditions
ux0 = @(t) pi*cos(pi*t);
ux1 = @(t) pi*cos(pi - pi*t);
u0 = @(x) sin(pi*x);
ut0 = @(x) -pi*cos(pi*x);

% Discretization
h = 1/N;                  % Element size
x = (0:h:1)';             % Node coordinates
numNodes = N + 1;         % Number of nodes
numTimeSteps = round(T/dt);

% Assemble the mass matrix and stiffness matrix
M = MassMatrix(N+1);
A = StiffnessMatrix(N, c);

% Initialize solution vectors
u = u0(x);               % Initial condition
ut = ut0(x);             % Initial time derivative

% Time integration using second-order Runge-Kutta method
for n = 1:numTimeSteps
    t = n * dt;              % Current time
    
    % Runge-Kutta step 1
    k1 = dt * (A*u - (M*(ut/dt)));
    u_mid = u + 0.5 * k1;
    
    % Runge-Kutta step 2
    k2 = dt * (A*u_mid - (M*(ut/dt)));
    u = u + k2;
    
    % Update time derivative
    ut = (u - u_mid) / dt;
end

% Plot the approximate and exact solutions
plot(x, u, 'b-', x, u_exact, 'r--');
legend('Approximate solution', 'Exact solution');
xlabel('x');
ylabel('u');
title(['Solution at t = ', num2str(T)]);


