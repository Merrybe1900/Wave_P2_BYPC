function residual = waveBackwardSolverP2massLumping2(meshT,bdaT, uhS,uDir,...
                                     u,u0,A,As,M,Me,ML,Nt,dt,displayOption)
%Note
% u0, to u3 corresponds to solution at Nt,Nt-1,...,Nt-3
residual = zeros(meshT.Ndof,1);
Ad = A(bdaT.freenode,bdaT.freenode);
Md = M(bdaT.freenode,bdaT.freenode);
MLd = spdiags(Md,0);
MLd = 1./MLd;

% ending zeros conditions
z0 = zeros(meshT.Ndof,1);
z1 = zeros(meshT.Ndof,1);



%%  Backward solver %
lambda = [1,0,0.5;0,1,0.5]';
weight = [1/6,1/6,2/3]; 

psi = zeros(3,3);
psi(:,1) = (2*lambda(:,1)-1).*lambda(:,1);
psi(:,2) = (2*lambda(:,2)-1).*lambda(:,2);
psi(:,3) = 4*lambda(:,1).*lambda(:,2);
bIn = zeros(length(bdaT.EdgeSOR)-1,3);

for n = Nt-2:-1:2 % first order ABC to initialize
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
     
    bz = accumarray([bdaT.EdgeSOR(1:end-1);bdaT.EdgeSOR(2:end);bdaT.SOR2node],bIn(:),[meshT.Ndof 1]); 
    z2 =  (bz  + (2*M-A)*z1 + (0.5*Me - M)*z0).*ML;
      
    u1 = zeros(meshT.Ndof,1);
    u1(bdaT.Dir) = uDir(:,n);
    ut = zeros(meshT.Ndof,1);
    ut(bdaT.Dir) = u(bdaT.Dir);
    b = (-A + 2*M)*ut;
    ut(bdaT.Dir) = u0(bdaT.Dir);
    b = b - M*(ut+u1); 
    b(bdaT.Dir) = u(bdaT.Dir);
     
    u1(bdaT.freenode) = (b(bdaT.freenode)-Ad*u(bdaT.freenode)...
                      + Md*(2*u(bdaT.freenode)-u0(bdaT.freenode)-u1(bdaT.freenode))).*MLd; 
    
    residual = residual +(u1-2*u+u0).*z2;
    
    if displayOption.displayB
        subplot(1,2,1)
        showsolution(meshT.node,meshT.elem,z2(1:meshT.nC))
        axis(displayOption.axis)
        view(displayOption.viewAngle)
        subplot(1,2,2)
        showsolution(meshT.node,meshT.elem,u1(1:meshT.nC))
        view(displayOption.viewAngle)
        axis(displayOption.axis)
        pause(0.01)
    end
    
    z0 = z1;
    z1 = z2;
    
    % update primal solution
    u0 = u;
    u = u1;
 end

residual = - residual/dt;
end
    


    
    
                          
                