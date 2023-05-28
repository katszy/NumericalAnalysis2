% Define the problem parameters
N = 4;
c = @ (x, t) 1 ;
f = @ (x, t) 1;
BC_Types = ["Neumann", "Neumann"];
BC_Values = {@(t) 0, @(t) 0};
T_final = 1;
num_time_steps = 100;

%
%%
c = @ (x , t ) 1;
f = @ (x , t ) 0;
g = @(x) sin(pi*x);
N = 4;
BC_Types = ["Dirichlet", "Dirichlet"];
BC_Values = {@(t) 0, @(t) 0};

%%
theta = 0;
[T, U] =FEM_HeatEq_1D(N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps, theta);
fprintf('theta = 0')
X = linspace (0 , 1 , N + 1);
mesh (X , T , U );
%%
theta=0.5;
[T, U] = FEM_HeatEq_1D(N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps, theta);
fprintf('theta = 0.5')
X = linspace (0 , 1 , N + 1);
mesh (X , T , U );

%%
theta=1;
[T, U] = FEM_HeatEq_1D(N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps, theta);
fprintf('theta =1')
X = linspace (0 , 1 , N + 1);
mesh (X , T , U );
%%
X = linspace (0 , 1 , N + 1);
mesh (X , T , U );

%%
% Compute and print the convergence table
u_exact = @(x, t) exp(-pi^2*t)*sin(pi*x);
L = [2, 6];
[e1, r1] = PrintConvergenceTableTheta(u_exact, L, 0, N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps);
[e2, r2] = PrintConvergenceTableTheta(u_exact, L, 0.5, N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps);
[e3, r3] = PrintConvergenceTableTheta(u_exact, L, 1, N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps);


%%
theta = 0.5 ;
f = @ (x , t ) 3 * x * exp ( -3* t ); 
c = @ (x , t ) 1;
N = 4;
BC_Types(1) = "Dirichlet";
BC_Types(2) = "Neumann";
BC_Values{2} = @ ( t ) (1 - exp ( -3* t ));
[T, U] = FEM_HeatEq_1D(N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps, theta);
X = linspace (0 , 1 , N + 1);
mesh (X , T , U );
%%
% Define the domain
x_min = 0;
x_max = 1;
t_min = 0;
t_max = 1;

% Define the grid resolution
Nx = 100; % Number of points in x-direction
Nt = 100; % Number of points in t-direction

% Generate the grid of points
[x, t] = meshgrid(linspace(x_min, x_max, Nx), linspace(t_min, t_max, Nt));

% Define the exact solution function
u_exact = @(x, t) x .* (1 - exp(-3*t));

% Evaluate the exact solution at the grid points
u = u_exact(x, t);

% Plot the surface
figure;
surf(x, t, u);
title('Exact Solution');
xlabel('x');
ylabel('t');
zlabel('u');
