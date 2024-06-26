 function [ e, r ] = PrintConvergenceTableTheta(u, L, theta, N, c, f, g, BC_Types, BC_Values, T_final, num_time_steps)
        %PRINTCONVERGENCETABLE This funcion computes and tabulates FEM errors
        %   Computes estimated errors for 1/k = 4, 8, 16, 32, and 64
        s = L(2) - L(1);
        e = zeros(s, 1);
        for k = 1:s+1
            l = L(1) + k;
            num_steps = 2^l;
            % Solve and compute error
            [~, Uh] = FEM_HeatEq_1D(num_steps, c, f, g, BC_Types, BC_Values,  T_final, num_time_steps, theta);
            uh = Uh(end, :)';
            e(k) = ErrorNorm(u, uh);
        end
        r = log2(e(1:(end-1)) ./ e(2:end));
        fprintf("\n");
        fprintf(join(repmat("%12s|", 1,s+1), ''), ["1/k" 2.^(1+(L(1):L(2)))]);
        fprintf("\n%12s|" + join(repmat("%12.4e|",1,s+1), ""), ["error"; e]);
        fprintf("\n%12s|" + join(repmat("%12.4f|",1,s+1), ""), ["rate"; "n/a"; r]);
        fprintf("\n");
 end 