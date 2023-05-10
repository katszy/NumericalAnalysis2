function [result, interval] = adint(f, left, right, M, p, w, tol, counter)
%performs adaptive numerical integration on interval [a, b] on function f
%M: number of quadrature points 
%p: order
%w: raw weights
%tol: tolerance 

global interval;
if counter==1
    interval=1;
else
    interval = interval+1;
end
result=0;
mid = (right+left)/2;

%calulate Q
step = (right-left)/(M-1);
Q=0;
xi=left;
for j=1:M
    Q = Q + f(xi)*w(j);
    xi=xi+step;
end
Q = Q*step;

%calculate Q1
step = (mid-left)/(M-1);
Q1=0;
xi=left;
for j=1:M
    Q1 = Q1 + f(xi)*w(j);
    xi=xi+step;
end
Q1 = Q1*step;

%calculate Q2
step = (right-mid)/(M-1);
Q2=0;
xi=mid;
for j=1:M
    xi=xi+step;
    Q2 = Q2 + f(xi)*w(j);
end
Q2 = Q2*step;

h = abs((Q1+Q2-Q)/(2.^p-1));

if h < tol
    result = Q1+Q2;
else
    result = adint(f, left, (left+right)/2, M, p, w, tol/2, counter+1) + adint(f, (left+right)/2, right, M, p, w, tol/2, counter+1);
end

if counter == 1
    result = result*(length(w)-1);
end

return

end

