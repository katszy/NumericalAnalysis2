function F = lotka_volterra(a, b, c, d, v)
    F(1) = a*v(1)-b*v(1)*v(2); 
    F(2) = c*v(1)*v(2)-d*v(2);
end

