function [] = solve(N,c,f,u,left_type, left_value, right_type, right_value)

% construct stiffness matrix, compute load vector, apply boundary conditions
[A, b] = apply_boundary_conditions(stiffnessmatrix(N, c),load_vector(N, f), left_type, left_value, right_type, right_value);
xi = A\b;

% plot solutions  
x1 = 0:0.01:1;
x2 = linspace(0, 1, N);
figure;
hold on;
plot(x1, u(x1)); % exact solution 
plot(x2, xi); % numerical solution
title("Solution of the Finite Element Method");
xlabel("x");
ylabel("y");
legend('Exact solution', 'Numerical solution', 'Location','northwest');

% calculate error
intervals = [4, 8, 16, 32, 64];
error = zeros(1, length(intervals));
for i = 1:length(intervals)
    N = intervals(i);
    [A, b] = apply_boundary_conditions(stiffnessmatrix(N, c),load_vector(N, f), left_type, left_value, right_type, right_value);
    xi = A\b;
    error(i) = errornorm(u,xi);
end

% calculate convergence rate -> second order
convergence_rate = zeros(1, length(intervals));
convergence_rate(1) = NaN;
for i=2:length(intervals)
    convergence_rate(i) = error(i-1)/error(i);
end

% plot 
figure;
hold on;
plot(intervals, error);
title("Errors for N = 4, 8, 16, 32, 64")
xlabel("N")
ylabel("Error")


% display error and convergence rates 
fprintf('Intervals\tError\t\tConvergence Rate\n');
fprintf('-----------------------------------------------\n');
for i = 1:numel(intervals)
    fprintf('%d\t\t%s\t\t%s\n', intervals(i), num2str(error(i)), num2str(convergence_rate(i)));
end

