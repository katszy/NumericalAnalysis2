function [A, b] = apply_boundary_conditions(A, b, BC_Types, BC_Values,t)
if BC_Types(1) == "Neumann"
     b(1) = b(1) - BC_Values{1}(t);
elseif BC_Types(1) == "Dirichlet"
    A(1,1) = 1;
    A(1,2) = 0;
    b(1) =  BC_Values{1}(t);
end

if BC_Types(2) == "Neumann"
    b(end) = b(end) + BC_Values{2}(t);
elseif BC_Types(2) == "Dirichlet"
    A(end,end) = 1;
    A(end,end-1) = 0;
    b(end) = BC_Values{2}(t);
end

end

