function [A_matrix, b_vector] = apply_boundary_conditions(A_matrix, b_vector, left_type, left_value, right_type, right_value)
if left_type == "Neumann"
     b_vector(1) = b_vector(1) + -1*left_value;
elseif left_type == "Dirichlet"
    A_matrix(1,1) = 1;
    A_matrix(1,2) = 0;
    b_vector(1) = left_value;
end

if right_type == "Neumann"
    b_vector(end) = b_vector(end) + right_value;
elseif right_type == "Dirichlet"
    A_matrix(end,end) = 1;
    A_matrix(end,end-1) = 0;
    b_vector(end) = right_value;
end

end

