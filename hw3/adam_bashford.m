function res = adam_bashford(a, b, c, d, x0, y0, T, N)
    h=T/N; 
     
    res = zeros(2,N);   
    res(:,1) = [x0,y0];
   
    F_LV = @lotka_volterra;
    F_CN = @crank_nicholson;
    
    res(:, 2) = fsolve(@(x) F_CN(a, b, c, d, h, x, res(:,1)), [x0,y0]);
    
    for i=3:N
        t_n1 = F_LV(a, b, c, d, res(:,i-1));
        t_n  = F_LV(a, b, c, d, res(:,i-2));
        res(:, i) = res(:,i-1)+h*(3/2*t_n1-1/2*t_n)';
    end
         
end
