function [x,y] = eulerfunc(a, b, c, d, x0, y0, T, N)
  
    h = T / N;   
    x = zeros(N + 1, 1);
    y = zeros(N + 1, 1);

    x(1) = x0;
    y(1) = y0;

    % Euler method
    for i = 1:N
        x_prime = a * x(i) - b * x(i) * y(i);
        y_prime = c * x(i) * y(i) - d * y(i);
        x(i + 1) = x(i) + h * x_prime;
        y(i + 1) = y(i) + h * y_prime;
    end
   
end

