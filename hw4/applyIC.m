function [y] = applyIC(g,x)
 % Get initial values
    % This is done by interpolating the initial value function g
     n = length(x);
     y = zeros(n, 1);
    for i = 1:n
        y(i) = g(x(i));
    end
end

