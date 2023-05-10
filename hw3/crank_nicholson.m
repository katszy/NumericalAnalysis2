function F = crank_nicholson(a, b, c, d, h, v, u)
    F(1) = h/2*(a*v(1)-b*v(1)*v(2)+a*u(1)-b*u(1)*u(2))+u(1)-v(1);
    F(2) = h/2*(c*v(1)*v(2)-d*v(2)+c*u(1)*u(2)-d*u(2))+u(2)-v(2);
end