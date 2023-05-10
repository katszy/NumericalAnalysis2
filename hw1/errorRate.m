function [errors, h_vec] = errorRate(f, left, right, N, M, p, w)
%helper function for error rates

% vector for errors
errors = zeros(N,1);
%vector for h values 
h = zeros(N-1,1);

for i = 1:N
    h(i)=(2^i);
    result = numint(f, left, right, 2^i, M, p, w); %result of numerical integration 
    errors(i) = abs(29.858325395498671-result)/29.858325395498671;  % the error is (1-result) since the integrand is 1 
end

h_vec = h(1:N-1);
errors = errors(1:N-1);
end

