% exercise 1
N = 64;
c = @(x) 1;
f = @(x) 6*x;
u = @(x) -x.^3 + 3*x;

left_type = "Dirichlet";
left_value = 0;
right_type = "Neumann";
right_value = 0;

solve(N,c,f,u,left_type, left_value, right_type, right_value);

%%
% exercise 2
N = 64;
c = @(x) 1+x;
f = @(x) 4*x;
u = @(x) 2*x-x.^2+log(1+x)-log(2);

left_type = "Neumann";
left_value = 3;
right_type = "Dirichlet";
right_value = 1;

solve(N,c,f,u,left_type, left_value, right_type, right_value);