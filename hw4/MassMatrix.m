function M = MassMatrix(N)
    h = 1 / N;  
    main_diag = h * ones(N, 1);
    main_diag([1, N]) = h / 2;
    M = spdiags(main_diag, 0, N, N);
end
