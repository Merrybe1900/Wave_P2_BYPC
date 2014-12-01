function [A, As, area] = stiffMatrixWaveP2(meshT,bdaT,varargin)



%% Compute geometric quantities and gradient of local basis
[Dlambda,area] = gradbasis(meshT.node,meshT.elem);

%% Assemble stiffness matrix
% Order for mass lumping: 3.5
% Otherwise: 2
[lambda, w] = quadptsWave(meshT.quadorder);
nQuad = size(lambda,1);

if meshT.quadorder == 3.5
    sA = zeros(28*meshT.nE,nQuad);
    ii = zeros(28*meshT.nE,1); jj = zeros(28*meshT.nE,1); 
    index = 0;
    for i = 1:7
        for j = i:7
            ii(index+1:index+meshT.nE) = double(meshT.elem2dof(:,i)); 
            jj(index+1:index+meshT.nE) = double(meshT.elem2dof(:,j));  
            index = index + meshT.nE;
        end
    end
    
    for p = 1:nQuad
        Dphip(:,:,1) = (4*lambda(p,1)-1).*Dlambda(:,:,1) ...
                     + 3*(Dlambda(:,:,1)*lambda(p,2)*lambda(p,3) ...
                     + Dlambda(:,:,2)*lambda(p,1)*lambda(p,3)...
                     + Dlambda(:,:,3)*lambda(p,1)*lambda(p,2));  
        Dphip(:,:,2) = (4*lambda(p,2)-1).*Dlambda(:,:,2) ...
                     + 3*(Dlambda(:,:,1)*lambda(p,2)*lambda(p,3) ...
                     + Dlambda(:,:,2)*lambda(p,1)*lambda(p,3)...
                     + Dlambda(:,:,3)*lambda(p,1)*lambda(p,2));  
               
        Dphip(:,:,3) = (4*lambda(p,3)-1).*Dlambda(:,:,3) ...
                     + 3*(Dlambda(:,:,1)*lambda(p,2)*lambda(p,3) ...
                     + Dlambda(:,:,2)*lambda(p,1)*lambda(p,3)...
                     + Dlambda(:,:,3)*lambda(p,1)*lambda(p,2));  
         
        Dphip(:,:,4) = 4*(Dlambda(:,:,2)*lambda(p,3) + Dlambda(:,:,3)*lambda(p,2))...
                     - 12*(Dlambda(:,:,1)*lambda(p,2)*lambda(p,3)...
                     + Dlambda(:,:,2)*lambda(p,1)*lambda(p,3)... 
                     + Dlambda(:,:,3)*lambda(p,1)*lambda(p,2));
        Dphip(:,:,5) = 4*(Dlambda(:,:,1)*lambda(p,3) + Dlambda(:,:,3)*lambda(p,1))...
                     - 12*(Dlambda(:,:,1)*lambda(p,2)*lambda(p,3)...
                     + Dlambda(:,:,2)*lambda(p,1)*lambda(p,3)... 
                     + Dlambda(:,:,3)*lambda(p,1)*lambda(p,2));
        Dphip(:,:,6) = 4*(Dlambda(:,:,1)*lambda(p,2) + Dlambda(:,:,2)*lambda(p,1))...
                     - 12*(Dlambda(:,:,1)*lambda(p,2)*lambda(p,3)...
                     + Dlambda(:,:,2)*lambda(p,1)*lambda(p,3)... 
                     + Dlambda(:,:,3)*lambda(p,1)*lambda(p,2));                
        Dphip(:,:,7) = 27*(Dlambda(:,:,1)*lambda(p,2)*lambda(p,3)...
                     + Dlambda(:,:,2)*lambda(p,1)*lambda(p,3)... 
                     + Dlambda(:,:,3)*lambda(p,1)*lambda(p,2));
        index = 0;
        for i = 1:7
            for j = i:7
                Aij =  w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
                Aij = Aij.*area;
                sA(index+1:index+meshT.nE,p) = Aij;
                index = index + meshT.nE;
            end
        end
    end
    
    Kn = varargin{1};
    %As = stiffMatrixABCP2(meshT,bdaT,Kn);
    As =[];
else
    % compute non-zeros
    ii = zeros(21*meshT.nE,1); jj = zeros(21*meshT.nE,1); 
    index = 0;
    for i = 1:6
        for j = i:6
            ii(index+1:index+meshT.nE) = double(meshT.elem2dof(:,i)); 
            jj(index+1:index+meshT.nE) = double(meshT.elem2dof(:,j));  
            index = index + meshT.nE;
        end
    end
    sA = zeros(21*meshT.nE,nQuad);
    for p = 1:nQuad
        % Dphi at quadrature points
        Dphip(:,:,1) = (4*lambda(p,1)-1).*Dlambda(:,:,1);            
        Dphip(:,:,2) = (4*lambda(p,2)-1).*Dlambda(:,:,2);            
        Dphip(:,:,3) = (4*lambda(p,3)-1).*Dlambda(:,:,3);            
        Dphip(:,:,4) = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
        Dphip(:,:,5) = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
        Dphip(:,:,6) = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
        index = 0;
        for i = 1:6
            for j = i:6

                Aij =  w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
                Aij = Aij.*area;
                sA(index+1:index+meshT.nE,p) = Aij;
                index = index + meshT.nE;
            end
        end
    end
    As = stiffMatrixABCP2(meshT,bdaT);
end
sA = sum(sA,2);
% assemble the matrix
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),meshT.Ndof,meshT.Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),meshT.Ndof,meshT.Ndof);
A = A + AU + AU';
clear Aij ii jj sA

end



