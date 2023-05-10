function [b] = load_vector(N,f)
h = 1/N;
b = zeros(N,1);
for i = 1:(N-1)
    b(i) = h*f((i*h));
end
b(1) = (h/2)*f(0);
b(N) = h/2*f(1);
end
  