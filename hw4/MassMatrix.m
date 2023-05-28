function M = MassMatrix(N)
    h = 1 / N;  % Subinterval width
    main_diag = h * ones(N, 1);
    
    % Update the diagonal entries
    main_diag([1, N]) = h / 2;
    
    % Construct the mass matrix using spdiags
    M = spdiags(main_diag, 0, N, N);
end
