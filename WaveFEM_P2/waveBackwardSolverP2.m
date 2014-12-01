function residual = waveBackwardSolverP2(meshT,bdaT,uhS,uDir,...
                      u,u0,A,As,M,Me,Ms,Nt,dt,displayOption)
%***************** Backward wave classcial function *****************%
% Author: Yun Bai
% Date: 21.08.2014
%********************************************************************%
% Note: For the time being, the auxiliary absorbing boundary function is 
%       not working properly, therefore disabled.
%********************************************************************%

% Ending zero conditions for dual problem %
z0 = zeros(meshT.Ndof,1);
z = zeros(meshT.Ndof,1);


% Initial absorbing boundary condition
% phi = cell(meshT.kABC,1);
% phid = phi;
% phidd = phi;
% for i = 1:meshT.kABC
%     phi{i} = zeros(2*(2*bdaT.NS(i)-1),1);
%     phid{i} = zeros(2*(2*bdaT.NS(i)-1),1);
%     phidd{i} = zeros(2*(2*bdaT.NS(i)-1),1);
% end
% 
% zs = cell(meshT.kABC,1);
% zt = cell(meshT.kABC,1);
% Kt = cell(meshT.kABC,1);

%--- Backward solver ---%
residual = zeros(meshT.Ndof,1);
Ad = A(bdaT.freenode,bdaT.freenode);
Md = M(bdaT.freenode,bdaT.freenode);


[lambda,weight] = quadpts1(4);
psi = zeros(3,3);
psi(:,1) = (2*lambda(:,1)-1).*lambda(:,1);
psi(:,2) = (2*lambda(:,2)-1).*lambda(:,2);
psi(:,3) = 4*lambda(:,1).*lambda(:,2);

bIn = zeros(bdaT.lenSOR-1,3);

for n=Nt:-1:3 
     
    %---- Dual problem ----%
    zIn1 = uhS(1:bdaT.lenSOR-1,n);
    zIn2 = uhS(2:bdaT.lenSOR,n);
    zIn3 = uhS(bdaT.lenSOR+1: 2*bdaT.lenSOR-1,n);
    
    % Neumann source boundary condition
    bIn(:,1) = (weight(1)*psi(1,1)*zIn1 + weight(2)*psi(2,1)*zIn2...
             + weight(3)*psi(3,1)*zIn3).*bdaT.elSOR;
    bIn(:,2) = (weight(1)*psi(1,2)*zIn1 + weight(2)*psi(2,2)*zIn2...
             + weight(3)*psi(3,2)*zIn3).*bdaT.elSOR;
    bIn(:,3) = (weight(1)*psi(1,3)*zIn1 + weight(2)*psi(2,3)*zIn2...
             + weight(3)*psi(3,3)*zIn3).*bdaT.elSOR;
     
     b = accumarray([bdaT.EdgeSOR(1:end-1);bdaT.EdgeSOR(2:end);bdaT.SOR2node],bIn(:),[meshT.Ndof 1]); 
    % absorbing boundary condition
%     for i = 1:meshT.kABC
%         zs{i} = [z0(bdaT.EdgeABC{i});z0(bdaT.ABCid{i}+ meshT.nC)];
%         zt{i} =[(z(bdaT.EdgeABC{i}(1))- z0(bdaT.EdgeABC{i}(1)))/dt,...
%                 (z(bdaT.EdgeABC{i}(bdaT.NS(i)))- z0(bdaT.EdgeABC{i}(bdaT.NS(i))))/dt];
%         Kt{i} = [sqrt(Kn(bdaT.EdgeABC{i}(1))), sqrt(Kn(bdaT.EdgeABC{i}(bdaT.NS(i))))];
%     end
% 
%     [phi,phid,phidd,Aux] =  WaveABCP2(zs,zt,phi,phid,phidd,Ms,As,...
%                                       bdaT.NS,bdaT.bd2node,Kt,...
%                                       -dt,meshT.Ndof,meshT.kABC);
% 
%     z1 = (M+0.5*Me)\((b-A*z) + M *(2*z-z0) + 0.5*Me*z0 - Me*Aux*dt);
    
    % solving the dual problem
    z1 = (M+0.5*Me)\((b-A*z) + M *(2*z-z0) + 0.5*Me*z0);
    z0 = z;
    z = z1;
    

    %---- Recover Primal solution ----%
    % Dirichlet boundary condition

    u1 = zeros(meshT.Ndof,1);
    u1(bdaT.Dir) = uDir(:,n-2);
    ut = zeros(meshT.Ndof,1);
    ut(bdaT.Dir) = u(bdaT.Dir);
    b = -A*ut + 2*M*ut;
    ut(bdaT.Dir) = u0(bdaT.Dir);
    b = b -M*ut-M*u1; 
    b(bdaT.Dir) = u(bdaT.Dir);
     
    % solve the primal problem    
    u1(bdaT.freenode) = Md\(b(bdaT.freenode)-Ad*u(bdaT.freenode)...
                 + Md*(2*u(bdaT.freenode)-u0(bdaT.freenode))...
                 - Md*u1(bdaT.freenode)); 
       
    %---Estimate the residual at time n---% 
    residual = residual +(u1-2*u+u0).*z/dt;
     
    % update primal solution
    u0 = u;
    u = u1;
    
    if displayOption.displayB
        subplot(1,2,1)
        showsolution(meshT.node,meshT.elem,z(1:meshT.nC))
        axis(displayOption.axis)
        view(displayOption.viewAngle)
        subplot(1,2,2)
        showsolution(meshT.node,meshT.elem,u(1:meshT.nC))
        view(displayOption.viewAngle)
        axis(displayOption.axis)
        pause(0.01)
    end
end

residual = - residual;
end
    


    
    
        