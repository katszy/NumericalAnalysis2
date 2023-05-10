function e = errornorm(u_e,u_h)
% ERRORNORM This function computes L2 error
% Computes an estimated distance between the
% exact solution and the numerical solution .
% The integral is approximated with three - point
% Gauss - Legendre quadrature rule
% u_e : Exact solution ( function )
% u_h : Numerical solution ( vector of nodal values )
n = length ( u_h ); % number of subintervals
h = 1.0 / ( n - 1); % length of subintervals
% rule for reference interval [0 , 1]
t = [ (1 - sqrt (3/5)) / 2 ; 1 / 2 ; (1 + sqrt (3/5)) / 2 ];
w = [ 5/18 ; 8/18 ; 5/18 ];
e2 = 0;
for i = 1:( n -1)
% evaluate at quadrature points
v_e = u_e ( (i -1) * h + h * t );
v_h = (1 - t ) * u_h ( i ) + t * u_h ( i +1);
% square difference and sum with weights
e2 = e2 + sum ( ( v_e - v_h ).^2 .* w );
end
e = sqrt ( h * abs( e2 ) );
end

