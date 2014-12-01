function As = stiffMatrixABCP2(meshT,bdaT,varargin)

order = 2;

if meshT.quadorder ==3.5 % mass lumping ABC boundary condition (global boundary stiff matrix)
    
   Kn = varargin{1};
   
   elem2dof =[];
   Ka = [];
   for k = 1:meshT.kABC
       Kt = 1./sqrt(Kn(bdaT.EdgeABC{k}));
       elem2doft = [bdaT.EdgeABC{k}(1:end-1),bdaT.EdgeABC{k}(2:end), bdaT.ABCid{k}+ meshT.nC];
       elem2dof = [elem2dof;elem2doft];
       Ka = [Ka; 0.25*(Kt(1:end-1) +Kt(2:end)).*bdaT.elABC{k}];
   end
  
   elem2dof = double(elem2dof);
   lambda = [1,0,0.5; 0,1,0.5]';
   weight = [1/6,1/6,2/3];
   nQuad = size(lambda,1);
   elABC = cell2mat(bdaT.elABC);
   Dlambda = [-1./elABC,1./elABC];
   As = sparse(meshT.Ndof,meshT.Ndof);
   
   for p = 1:nQuad
        Dphip(:,1) = (4*lambda(p,1)-1).*Dlambda(:,1);
        Dphip(:,2) = (4*lambda(p,2)-1).*Dlambda(:,2);
        Dphip(:,3) = 4*(lambda(p,1)*Dlambda(:,2) + lambda(p,2)*Dlambda(:,1));

        for i = 1:3
            for j = 1:3   % only compute half of the off-diagonal part (local boundary stiff matrix)
                Aij = weight(p)*Dphip(:,i).*Dphip(:,j).*Ka;
                As = As + sparse(elem2dof(:,i),elem2dof(:,j),Aij,meshT.Ndof,meshT.Ndof);
            end
        end
       
   end
   
else  % classical boundary stiffness matrix for ABC

    As = cell(meshT.kABC,1);
    elABC = bdaT.elABC;
    
    parfor k = 1: meshT.kABC
            nbd = length(elABC{k});
            elem2dof = [(1:nbd)',(2:nbd+1)',((nbd+2):(2*nbd+1))'];
            Ndof = 2*nbd+1;

            % quadratic bases (1---3---2)
            [lambda,weight] = quadpts1(order);
            nQuad = size(lambda,1);
            Dlambda = [-1./elABC{k},1./elABC{k}];

            Dphip = zeros(nbd,3);
            A = sparse(Ndof,Ndof);
            for p = 1:nQuad
                Dphip(:,1) = (4*lambda(p,1)-1).*Dlambda(:,1);
                Dphip(:,2) = (4*lambda(p,2)-1).*Dlambda(:,2);
                Dphip(:,3) = 4*(lambda(p,1)*Dlambda(:,2) + lambda(p,2)*Dlambda(:,1));

                for i = 1:3
                    for j = 1:3   % only compute half of the off-diagonal part
                        Aij = weight(p)*Dphip(:,i).*Dphip(:,j).*elABC{k};
                        A = A + sparse(elem2dof(:,i),elem2dof(:,j),Aij,Ndof,Ndof);
                    end
                end
            end

            As{k} = A;
    end
end

end