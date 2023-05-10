f = @(x) 1./(0.01 + (x - 0.3).^2) + 1./(0.04 + (x - 0.9).^2) - 6;
a = 0;
b = 1;
N = 1000;
ground_truth  = 29.858325395498671;
fprintf('Ground truth 29.858325395498671 \n \n');

% equidistant numerical integration 
result1 = numint(f, a, b, N, 2, 2, [1/2,1/2]);
result2 = numint(f, a, b, N, 3, 4, [1/6,2/3,1/6]);
result3 = numint(f, a, b, N, 4, 4, [1/8,3/8,3/8,1/8]);
result4 = numint(f, a, b, N, 5, 6, [7/90,16/45,4/30,16/45,7/90]);

fprintf('Equidistant integration for rule 1: %.6f\n', result1);
fprintf('Equidistant integration for rule 2: %.6f\n', result2);
fprintf('Equidistant integration for rule 3: %.6f\n', result3);
fprintf('Equidistant integration for rule 4: %.6f\n \n', result4);

%adaptive numerical integration
tol=10e-4;
counter  = 1;
[result11, interval1] = adint(f, a, b, 2, 2, [1/2,1/2], tol, counter);
[result22, interval2] = adint(f, a, b, 3, 4, [1/6,2/3,1/6], tol, counter);
[result33, interval3] = adint(f, a, b, 4, 4, [1/8,3/8,3/8,1/8], tol, counter);
[result44, interval4] = adint(f, a, b, 5, 6, [7/90,16/45,4/30,16/45,7/90], tol, counter);
fprintf('Adaptive integration for rule 1: %.6f, intervals: %.0f \n', result11, interval1);
fprintf('Adaptive integration for rule 2: %.6f, intervals: %.0f\n', result22,interval2);
fprintf('Adaptive integration for rule 3: %.6f, intervals: %.0f\n', result33,interval3);
fprintf('Adaptive integration for rule 4: %.6f, intervals: %.0f\n \n', result44,interval4);


%let's use  M=2
tol = [];
for i=0:4
    tol = 1/10^i;
    counter = 1;
    [res, interval] = adint(f, a, b, 2, 2, [1/2,1/2], tol, counter);
    fprintf('eps = %e, #intervals.= %.0f, error rate %.6f \n', tol, interval, abs(res-29.858325395498671)/29.858325395498671)
end
% okay 10e-4 is good enough 
% now let's find an N value that produces a similar result
I = 20; %well actually 2^i = 1..N 
[errors, h] = errorRate(f, a, b, I, 2, 2, [1/2,1/2]); 
for i=1:length(errors)-1
    if errors(i) < 0.000001 
        fprintf('N: %.0f, error rate: %.6f \n',h(i), errors(i));
        break
    end
end

%comparing running times for M=2
fprintf('\n Running time for equidistant method, M=2: \n');
tic
result2 = numint(f, a, b, N, 3, 4, [1/6,2/3,1/6]);
toc

fprintf(' \n Running time for adaptive method, M=2: \n');
tic
[result22, interval2] = adint(f, a, b, 3, 4, [1/6,2/3,1/6], tol, counter); % oof very slow :(
toc
 

%2 point Gaussian quadrature 
N=100;
intervals = linspace(0,1,N+1);
res = 0;
for i=1:N
    res = res + gaussian2(f,intervals(i),intervals(i+1));
end
fprintf('\n');
fprintf('2 point Gaussian quadrature rule: %.6f \n', res);

for N=1:6
    intervals = linspace(0,1,2^N+1);
    res = 0;
    tic
    for i=1:length(intervals)-1
        res = res + gaussian2(f,intervals(i),intervals(i+1));
    end
    fprintf('\n N = %.0f, error rate %.8f \n', 2^N, abs(res-29.858325395498671)/29.858325395498671);
    toc 
end

fprintf('Similar error rates eg. p=4, p=6, but with the right parameters all of them produced pretty similar ones');
