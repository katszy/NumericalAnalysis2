function [] = convergence_rate(ref, a, b, c, d, x0, y0, T, method)
    error = zeros(1, 8);
    
    if method=="Euler"
        for i = 10:17
            N = 2^i;
            [res_x, res_y] = eulerfunc(a, b, c, d, x0, y0, T, N);
            calc = res_x(end) + res_y(end);
            error(i-9) = abs(ref - calc) / ref;
        end
    end
    
    if method=="AB"
        for i = 10:17
        	N = 2^i;
            res = adam_bashford(a, b, c, d, x0, y0, T, N);
            calc = res(1,end) + res(2,end);
            error(i-9) = abs(ref - calc) / ref;
        end
    end
      
%     % calculate convergence rate
%     convergence_rate = zeros(1, 8);
%     for i = 2:8
%         convergence_rate(i) = (error(i-1) / error(i));
%     end

    figure;
    plot(1:8, error(1:8));
    xlabel('Iteration');
    ylabel('Convergence rate');
    title(['Convergence rate vs. iteration for method ', method]);
end

