function [M,Me,varargout] = massMatrixWaveP2(meshT,bdaT,area,Ke,Kn,dt)

%********** massMatrixWaveP2 *******************%
% This function assembles mass matrix for either
% p2 masslumping domain mass matrix and global boundary mass matrix
% or classical p2 elements domain mass matrix, global boundary mass matrix
% and local boundary mass matrix.
%***********************************************%


% mass matrix in the computation domain

if meshT.quadorder == 3.5 % mass lumping mass matrix
    
    Ka = area.*Ke/dt/dt;
    Md = [1/20*Ka,1/20*Ka,1/20*Ka,...
          2/15*Ka,2/15*Ka,2/15*Ka,9/20*Ka];
    M = accumarray(meshT.elem2dof(:),Md(:),[meshT.Ndof 1]); 
    Me =  massMatrixBd(meshT,bdaT,Kn,dt);
    ML = M + 0.5*spdiags(Me,0);
    M = spdiags(M,0,meshT.Ndof,meshT.Ndof);
    varargout{1} = 1./ML;
    
else % classical mass matrix
    
    ii = zeros(21*meshT.nE,1); jj = zeros(21*meshT.nE,1);
    index = 0;
    
    for i = 1:6
        for j = i:6
            ii(index+1:index + meshT.nE) = double(meshT.elem2dof(:,i)); 
            jj(index+1:index + meshT.nE) = double(meshT.elem2dof(:,j));  
            index = index + meshT.nE;
        end
    end
    
    [lambda,weight] = quadpts(meshT.quadorder+2);
    phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
    phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
    phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
    phi(:,4) = 4*lambda(:,2).*lambda(:,3);
    phi(:,5) = 4*lambda(:,3).*lambda(:,1);
    phi(:,6) = 4*lambda(:,1).*lambda(:,2);
   
    sA = zeros(21*meshT.nE,1);
    index = 0;
    for i = 1:6
         for j = i:6
             aij =  Ke*(weight*(phi(:,i).*phi(:,j))); % use piecewise constant...
                                                    % to simulate coefficient
             Aij = aij.*area/dt/dt;
             sA(index+1:index + meshT.nE) = Aij;
             index = index + meshT.nE;
         end
    end

    % assemble the matrix
    diagIdx = (ii == jj);   upperIdx = ~diagIdx;
    M = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),meshT.Ndof,meshT.Ndof);
    AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),meshT.Ndof,meshT.Ndof);
    M = M + AU + AU';
    [Me, Ms] =  massMatrixBd(meshT,bdaT,Kn,dt);
    varargout{1} = Ms;
end
end


%%
function [Me, varargout] =  massMatrixBd(meshT,bdaT,Kn,dt)
% assemble the boundary mass matrix
elem2dof = [];
for k = 1 : meshT.kABC
    elem2doft = [bdaT.EdgeABC{k}(1:end-1), bdaT.EdgeABC{k}(2:end), bdaT.ABCid{k}+ meshT.nC];
    elem2dof = [elem2dof;elem2doft];
end
elem2dof = double(elem2dof);
elABCt = cell2mat(bdaT.elABC);
if meshT.quadorder == 3.5
    lambda = [1,0,0.5;0,1,0.5]';
    weight = [1/6,1/6,2/3];
else
    order = 4;
    [lambda,weight] = quadpts1(order);
    weight = weight';
end
nQuad = size(lambda,1);
Me = sparse(meshT.Ndof,meshT.Ndof);
phi = zeros(nQuad,3);
phi(:,1) = 2*lambda(:,1).^2 - lambda(:,1);
phi(:,2) = 2*lambda(:,2).^2 - lambda(:,2);
phi(:,3) = 4*lambda(:,1).*lambda(:,2);

Kt = (sqrt(Kn(elem2dof(:,1))) +  sqrt(Kn(elem2dof(:,2))))*0.5;
for i = 1:3
    for j = 1:3   % only compute half of the off-diagonal part (local boundary stiff matrix)
        Mij = weight*(phi(:,i).*phi(:,j))*elABCt.*Kt/dt;
        Me = Me + sparse(elem2dof(:,i),elem2dof(:,j),Mij,meshT.Ndof,meshT.Ndof);
    end
end


if meshT.quadorder ~= 3.5
    Ms = cell(meshT.kABC,1);
    Kt = cell(meshT.kABC,1);
    
    for k = 1 : meshT.kABC
        Kt{k} = (Kn(bdaT.EdgeABC{k}(1:end-1)) + Kn(bdaT.EdgeABC{k}(2:end)))*0.5; 
    end
    
    elABC = bdaT.elABC;
    
    for k = 1 : meshT.kABC
        nbd = length(elABC{k});
        elem2dof = [(1:nbd)',(2:nbd+1)',((nbd+2):(2*nbd+1))'];
        elem2dof = double(elem2dof);
        Ndof = 2*nbd+1;

        % quadratic bases (1---3---2)
        [lambda,weight] = quadpts1(order);
        phi = zeros(nQuad,3);
        phi(:,1) = 2*lambda(:,1).^2 - lambda(:,1);
        phi(:,2) = 2*lambda(:,2).^2 - lambda(:,2);
        phi(:,3) = 4*lambda(:,1).*lambda(:,2);
        
        M = sparse(Ndof,Ndof);
        for i = 1:3
            for j = 1:3   % only compute half of the off-diagonal part (local boundary stiff matrix)
                Mij = weight'*(phi(:,i).*phi(:,j))*elABC{k}.*Kt{k};
                M = M + sparse(elem2dof(:,i),elem2dof(:,j),Mij,Ndof,Ndof);
            end
        end

        Ms{k} = M;
    end
    varargout{1} = Ms;
end

end


