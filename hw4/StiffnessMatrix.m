function A = StiffnessMatrix(N, c)
    h = 1 / N;
    A = zeros(N + 1, N + 1);
    A(1, 1) = (c(0) + c(h)) / (2 * h);
    A(1, 2) = -A(1, 1);
    for i = 2:N
        A(i, i - 1) = -(c(h * (i - 1)) + c(h * i)) / (2 * h);
        A(i, i) = (c(h * (i - 1)) + 2 * c(h * i) + c(h * (i + 1))) / (2 * h);
        A(i, i + 1) = -(c(h * i) + c(h * (i + 1))) / (2 * h);
    end
    A(N + 1, N) = -(c(1 - h) + c(1)) / (2 * h);
    A(N + 1, N + 1) = -A(N + 1, N);
end

