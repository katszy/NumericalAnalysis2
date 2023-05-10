function [result] = numint(f, left, right, N, M, p, w)
%performs equidistant numerical integration on interval [a, b] on function f
%M: number of quadrature points 
%p: order
%w: raw weights

h = (right - left) / N;
result=0;

for i=0:N-1
    a=left+i*h;
    b=a+h;
    
    step = (b-a)/M;
    xi=a;
    
    for j=1:M
        xi=xi+step;
        result = result + f(xi)*w(j);
    end
    
end

result=result*h;

end

