classdef HeatEquation
    %HEATEQUATION This class managed the data structures we use to solve the heat equation.
    
    properties
        % We always need to following paramets
        N   % number of intervals
        c   % diffusion coefficient (function) 
        f   % source term (function)
        
        % Default to homogeneous Neumann boundaries
        BC_Types  = ['N', 'N'];
        BC_Values = {@(t) 0, @(t) 0};
        
        % Default to zero initial condition
        g = @(x) 0*x;
    end
    
    methods
        % Constructor
        function obj = HeatEquation(N, c, f)
            obj.N = N;
            obj.c = c;
            obj.f = f;
        end
        
        function [ A ] = StiffnessMatrix(obj, t)
        %STIFNESSSSMATRIX Assemble stiffness matrix at time t
            A = StiffnessMatrix( obj.N, @(x) obj.c(x,t) );
        end

        function [ M ] = MassMatrix(obj)
        %MASSMATRIX Assemble mass matrix
            M = MassMatrix(obj.N);
        end
        
        function [ b ] = LoadVector(obj, t)
        %LOADVECTOR Assemble the load vector at time t
            b = LoadVector( obj.N, @(x) obj.f(x, t) );
            
            % Account for boundary conditions
            if obj.BC_Types(1) == 'N'
                b(1) = b(1) - obj.BC_Values{1}(t);
            end
            if obj.BC_Types(2) == 'N'
                b(end) = b(end) + obj.BC_Values{2}(t);
            end 
        end
        
        function [ A, b ] = ApplyDirichletBCs(obj, A, b, t)
        %APPLYDIRICHLETBCS Impose Dirichlet conditions on system
            if obj.BC_Types(1) == 'D'
                A(1,1) = 1; A(1,2) = 0;
                b(1) = obj.BC_Values{1}(t);
            end
            if obj.BC_Types(2) == 'D'
                A(end, end) = 1; A(end,end-1) = 0;
                b(end) = obj.BC_Values{2}(t);
            end
        end
        
        % Solve the heat equation with a theta scheme
        function [T, U] = SolveTheta(obj, T_final, num_time_steps, theta)
            % Assemble the mass matrix
            M = obj.MassMatrix();
            
            % Determine the length of the time steps
            n = obj.N + 1;                   % number of points in space
            k = T_final / num_time_steps;    % length of the time steps

            % Make a vector of the time partitioning
            T = linspace(0, T_final, num_time_steps + 1)';

            % Matrix to hold all the solutions
            U = zeros(num_time_steps + 1, n);
            
            % Get initial values
            % This is done by interpolating the initial value function g
            x = linspace(0, 1, n);
            u = zeros(n, 1);
            for i=1:n
                u(i) = obj.g(x(i));
            end
            U(1,:) = u;

            % Start time loop
            for l=1:num_time_steps
                % increase time to t_theta
                % FIXME: Put in the correct t_theta here
                t = k * l;
                
                % Get stiffness matrix A and loadvector b at time t=t_theta
                A = obj.StiffnessMatrix(t);
                b = obj.LoadVector(t);
                
                % Form the linear system to be solved
                AA = M + theta * k* A;
                bb = M * u + k * (1 - theta) * b;
                
                % Apply BCs (at time t_{l+1} = l * k)
                [AA, bb] = obj.ApplyDirichletBCs(AA, bb, l * k);
                
                % Solve for solution at time t and store it
                u = AA \ bb;
                U(l+1, :) = u';
            end         
        end
        
        %%function y = F(obj, t, u)
        %%F defines right hand side of ode system equation
        %    A = obj.StiffnessMatrix(t);
        %    b = obj.LoadVector(t);
        %
        %    y = 
        %end
        
        %%function [T, U] = SolveODE15s(obj, t_final)
        %%Assemble the mass matrix
        %M = obj.MassMatrix();
        %n = obj.N + 1;
        %
        %%Interpolate initial values
        %x = linspace(0, 1, n);
        %u0 = zeros(obj.N + 1, 1);
        %for i=1:n
        %    u0(i) = obj.g(x(i));
        %end
        %
        %%Solve
        %F = @obj.F;
        %JPattern = spdiags(ones(n, 3), [-1,0,1], n,n);
        %odeopts = 
        %[T, U] = 
        %end
        
        
        function [ e, r ] = PrintConvergenceTableTheta(obj, u, L, theta)
        %PRINTCONVERGENCETABLE This funcion computes and tabulates FEM errors
        %   Computes estimated errors for 1/k = 4, 8, 16, 32, and 64

        s = L(2) - L(1);
        e = zeros(s, 1);
        for k = 1:s+1
            l = L(1) + k;
            num_steps = 2^l;
            % Solve and compute error
            [~, Uh] = obj.SolveTheta(1, num_steps, theta); uh = Uh(end, :)';
            e(k) = Errornorm(u, uh);
        end
        r = log2(e(1:(end-1)) ./ e(2:end));
        fprintf(string('\n'));
        fprintf(join(repmat(string('%12s|'), 1,s+1), ''), [string('1/k') 2.^(1+(L(1):L(2)))]);
        fprintf(string('\n%12s|') + join(repmat(string('%12.4e|'),1,s+1), string('')), [string('error'); e]);
        fprintf(string('\n%12s|') + join(repmat(string('%12.4f|'),1,s+1), string('')), [string('rate'); string('n/a'); r]);
        fprintf(string('\n'));
        end
    end   
end