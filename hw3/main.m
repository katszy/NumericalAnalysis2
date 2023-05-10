a = 1;
b = 1;
c = 1;
d = 1;
x0 = 5;
y0 = 1;
T = 50;
h = 1e-2;
N = T/h;
%%
[x,y] = eulerfunc(a, b, c, d, x0, y0, T, N);

figure;
hold on;
plot(0:h:T, x');
plot(0:h:T, y');
legend('Rabbit', 'Fox');
xlabel('t');
ylabel('population');
title('Lotka-Volterra equation (Euler method)');
%%
ab_exact = adam_bashford(a, b, c, d, x0, y0, T, N);

x = ab_exact(1,:)';
y = ab_exact(2,:)';
figure;
hold on;
plot(0:h:T-h,x);
plot(0:h:T-h,y);
title('Lotka-Volterra equation (Adams-Bashford method)');
legend('Rabbit', 'Fox');
xlabel('t');
ylabel('population');

%%
N = 10^7;
T = 50;
ab_exact = adam_bashford(a, b, c, d, x0, y0, T, N); %exact solution for the a-b method, slow calculaion :( 
ref = ab_exact(1,end) + ab_exact(2,end);
[res_x, res_y] = eulerfunc(a, b, c, d, x0, y0, T, N); %exact solution for the euler method
%%
%interesting results
convergence_rate(ref, a, b, c, d, x0, y0, T, 'Euler'); 
convergence_rate(ref, a, b, c, d, x0, y0, T, 'AB');

%%
N = 10^7;
T = 50;
h = T/N;

fprintf("local maxima and minima values for Euler method \n");
exact_solution = res_x + res_y;
% plot(0:h:T,exact_solution');
% title('x(t)+y(t)');
find_peaks(exact_solution);


fprintf("local maxima and minima values for Adam-Bashford method \n");
%exact solution for the euler method is already calculated 
exact_solution_2 = ab_exact(1,:)' + ab_exact(2,:)';
find_peaks(exact_solution_2);
