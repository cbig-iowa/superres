function [D,Dt] = defDDt
        % defines finite difference operator D
        % and its transpose operator
        
        D = @(U) grad(U);
        Dt = @(V) div(V);
        
        function DU = grad(U)
            % Forward finite difference operator 
            %(with circular boundary conditions)
            DU(:,:,1) = [diff(U,1,2), U(:,1) - U(:,end)];
            DU(:,:,2) = [diff(U,1,1); U(1,:) - U(end,:)];
        end
        
        function DtXY = div(V)
            % Divergence operator (transpose of gradient)
            X = V(:,:,1);
            Y = V(:,:,2);
            DtXY = [X(:,end) - X(:,1), -diff(X,1,2)];
            DtXY = DtXY + [Y(end,:) - Y(1,:); -diff(Y,1,1)];
        end
end