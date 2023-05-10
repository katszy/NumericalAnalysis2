function I = gaussian2(f, a, b)
%2-point Gaussian quadrature rule
%f: function, [a,b] interval 

w = 1.0 * abs(b - a) / 2;
x = [-1/sqrt(3), 1/sqrt(3)];
%mapping [a, b] to [-1, 1]
t = (b - a)/2 .* x + (b + a)/2;
fx0=w*f(t(1));
fx1=w*f(t(2));
I=fx0+fx1;
end