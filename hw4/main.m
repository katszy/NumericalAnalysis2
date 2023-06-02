%% exact solutions  - plotted later 
x_min = 0;
x_max = 1;
t_min = 0;
t_max = 1; 
Nx = 100;
Nt = 100;

[x, t] = meshgrid(linspace(x_min, x_max, Nx), linspace(t_min, t_max, Nt));
u_exact1 = @(x, t) exp(-pi^2*t)*sin(pi*x);
u_exact2 = @(x, t) x .* (1 - exp(-3*t));

% Evaluate the exact solution at the grid points
u1 = u_exact1(x, t);
u2 = u_exact2(x, t);

%% EXERCISE 1
% Define the problem parameters
T_final = 1;
num_time_steps = 100;
c = @ (x , t ) 1;
f = @ (x , t ) 0;
g = @(x) sin(pi*x);
N = 4;
BC_Types = ["Dirichlet", "Dirichlet"];
BC_Values = {@(t) 0, @(t) 0};

figure;

% Compute and plot solution for theta = 0
subplot(2, 2, 1);
theta = 0;
[T, U] = FEM_HeatEq_1D(N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps, theta);
X = linspace(0, 1, N + 1);
mesh(X, T, U);
title('Theta = 0');
xlabel('x');
ylabel('t');
zlabel('u');

% Compute and plot solution for theta = 0.5
subplot(2, 2, 2);
theta = 0.5;
[T, U] = FEM_HeatEq_1D(N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps, theta);
X = linspace(0, 1, N + 1);
mesh(X, T, U);
title('Theta = 0.5');
xlabel('x');
ylabel('t');
zlabel('u');

% Compute and plot solution for theta = 1
subplot(2, 2, 3);
theta = 1;
[T, U] = FEM_HeatEq_1D(N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps, theta);
X = linspace(0, 1, N + 1);
mesh(X, T, U);
title('Theta = 1');
xlabel('x');
ylabel('t');
zlabel('u');

% Plot the exact solution
subplot(2, 2, 4);
mesh(x, t, u1);
title('Exact Solution');
xlabel('x');
ylabel('t');
zlabel('u');
hold on;

% Compute and print the convergence table
u_exact = @(x, t) exp(-pi^2*t)*sin(pi*x);
L = [2, 6];
[e1, r1] = PrintConvergenceTableTheta(u_exact, L, 0, N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps);
[e2, r2] = PrintConvergenceTableTheta(u_exact, L, 0.5, N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps);
[e3, r3] = PrintConvergenceTableTheta(u_exact, L, 1, N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps);
%% EXERCISE 2 

f = @(x, t) 3 * x * exp(-3 * t);
c = @(x, t) 1;
N = 4;
BC_Types(1) = "Dirichlet";
BC_Types(2) = "Neumann";
BC_Values{2} = @(t) (1 - exp(-3 * t));

figure;

% Compute and plot solution for theta = 0
subplot(2, 2, 1);
theta = 0;
[T, U] = FEM_HeatEq_1D(N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps, theta);
X = linspace(0, 1, N + 1);
mesh(X, T, U);
title('Theta = 0');
xlabel('x');
ylabel('t');
zlabel('u');

% Compute and plot solution for theta = 0.5
subplot(2, 2, 2);
theta = 0.5;
[T, U] = FEM_HeatEq_1D(N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps, theta);
X = linspace(0, 1, N + 1);
mesh(X, T, U); 
title('Theta = 0.5');
xlabel('x');
ylabel('t');
zlabel('u');

% Compute and plot solution for theta = 1
subplot(2, 2, 3);
theta = 1;
[T, U] = FEM_HeatEq_1D(N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps, theta);
X = linspace(0, 1, N + 1);
mesh(X, T, U);
title('Theta = 1');
xlabel('x');
ylabel('t');
zlabel('u');

% Plot the exact solution
subplot(2, 2, 4);
mesh(x, t, u2);
title('Exact Solution');
xlabel('x');
ylabel('t');
zlabel('u');


% Compute and print the convergence table
u_exact2 = @(x, t) x .* (1 - exp(-3*t));
L = [2, 6];
[e1, r1] = PrintConvergenceTableTheta(u_exact2, L, 0, N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps);
[e2, r2] = PrintConvergenceTableTheta(u_exact2, L, 0.5, N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps);
[e3, r3] = PrintConvergenceTableTheta(u_exact2, L, 1, N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps);

%% EXERCISE3 -> wave.m 