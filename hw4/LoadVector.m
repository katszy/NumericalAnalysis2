function [b] = LoadVector(N,f)
  h = 1 / N;
    b = zeros(N + 1, 1);
    b(1) = h / 2 * f(0);
    for i = 2:N
        b(i) = h * f(i * h);
    end
    b(N + 1) = h / 2 * f(1);
end
  