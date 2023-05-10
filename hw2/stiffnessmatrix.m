function A = stiffnessmatrix(N, c)
    h = 1 / (N-1);
    A = zeros(N, N);
    A(1, 1) = (c(0) + c(h)) / (2 * h);
    A(1, 2) = -(c(0) + c(h)) / (2 * h);
   for i = 2:(N-1)
       xi = i * h;
       xi_minus_1 = (i - 1) * h;
       xi_plus_1 = (i + 1) * h;
       A(i, i - 1) = -(c(xi_minus_1) + c(xi)) / (2 * h);                         
       A(i, i) = (c(xi_minus_1) + 2 * c(xi) + c(xi_plus_1)) / (2 * h);           
       A(i, i + 1) = -(c(xi) + c(xi_plus_1)) / (2 * h);
   end
    A(N, N-1) = -(c(1 - h) + c(1)) / (2 * h);
    A(N, N) = -A(N,N-1);
end

