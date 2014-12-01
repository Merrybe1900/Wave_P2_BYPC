function [uSOR, varargout] = waveForwardSolverP2massLumping(meshT,bdaT,source,signal,ks,A,As,M,Me,ML,Nt,dt,displayOption,Kn)
%**************** Forward wave mass lumping function ***************%
% Author: Yun Bai
% Date: 21.08.2014
%********************************************************************%

% Initial zeros conditions
u0 = zeros(meshT.Ndof,1);
u1 = zeros(meshT.Ndof,1);


%%  Forward solver %

lambda = [1,0,0.5;0,1,0.5]';
weight = [1/6,1/6,2/3];
pxy1 = meshT.node(bdaT.EdgeSOR(1:end-1),:);
pxy2 = meshT.node(bdaT.EdgeSOR(2:end),:);
pxy3 = 0.5*(pxy1 + pxy2);   

psi = zeros(3,3);
psi(:,1) = (2*lambda(:,1)-1).*lambda(:,1);
psi(:,2) = (2*lambda(:,2)-1).*lambda(:,2);
psi(:,3) = 4*lambda(:,1).*lambda(:,2);
bIn = zeros(length(bdaT.EdgeSOR)-1,3);
uSOR = zeros(length(bdaT.EdgeSOR)+length(bdaT.SOR2node),Nt);
pSORt = signal.pSOR(ks,:);

if meshT.solType == 'i'
   uDir = zeros(length(bdaT.Dir),Nt); 
end
    
tic
for n = 2:3 % first order ABC to initialize
    
    uIn1 = source(pxy1,pSORt,n*dt,signal);
    uIn2 = source(pxy2,pSORt,n*dt,signal);
    uIn3 = source(pxy3,pSORt,n*dt,signal);
          
  
    bIn(:,1) = (weight(1)*psi(1,1)*uIn1 + weight(2)*psi(2,1)*uIn2...
             + weight(3)*psi(3,1)*uIn3).*bdaT.elSOR;
    bIn(:,2) = (weight(1)*psi(1,2)*uIn1 + weight(2)*psi(2,2)*uIn2...
             + weight(3)*psi(3,2)*uIn3).*bdaT.elSOR;
    bIn(:,3) = (weight(1)*psi(1,3)*uIn1 + weight(2)*psi(2,3)*uIn2...
             + weight(3)*psi(3,3)*uIn3).*bdaT.elSOR;
         
    b = accumarray([bdaT.EdgeSOR(1:end-1);bdaT.EdgeSOR(2:end);bdaT.SOR2node],...
                    bIn(:),[meshT.Ndof 1]); 
     
    if n ==2
       u2 =  (b  + (2*M-A)*u1 + (0.5*Me - M)*u0).*ML;
       if meshT.solType == 'i'
          uDir(:,2) = u2(bdaT.Dir); 
       end
    elseif n == 3
       u3 =  (b  + (2*M-A)*u2 + (0.5*Me - M)*u1).*ML;
       if meshT.solType == 'i'
          uDir(:,3) = u3(bdaT.Dir); 
       end
    end
     
end

M1 = 2*M - A;
M2 = Me - 2*dt*As;
M3 = -0.5*Me + M;

for n=4:Nt-1  % secdon order mass lumption 

    uIn1 = (source(pxy1,pSORt,(n-1)*dt,signal)...
         -  source(pxy1,pSORt,(n-3)*dt,signal));
    uIn2 = (source(pxy2,pSORt,(n-1)*dt,signal)...
         -  source(pxy2,pSORt,(n-3)*dt,signal));
    uIn3 = (source(pxy3,pSORt,(n-1)*dt,signal)...
         -  source(pxy3,pSORt,(n-3)*dt,signal));

    bIn(:,1) = (weight(1)*psi(1,1)*uIn1 + weight(2)*psi(2,1)*uIn2...
             + weight(3)*psi(3,1)*uIn3).*bdaT.elSOR;
    bIn(:,2) = (weight(1)*psi(1,2)*uIn1 + weight(2)*psi(2,2)*uIn2...
             + weight(3)*psi(3,2)*uIn3).*bdaT.elSOR;
    bIn(:,3) = (weight(1)*psi(1,3)*uIn1 + weight(2)*psi(2,3)*uIn2...
             + weight(3)*psi(3,3)*uIn3).*bdaT.elSOR;


    b = accumarray([bdaT.EdgeSOR(1:end-1);bdaT.EdgeSOR(2:end);bdaT.SOR2node],...
                    bIn(:),[meshT.Ndof 1]); 
    % wave forward solver
    %u4 =  (b  + (2*M-A)*u3 + (Me - 2*dt*As)*u2 + (A-2*M)*u1 - (0.5*Me-M)*u0).*ML;
    u4 = (b + M1*(u3-u1) + M2*u2 + M3*u0).*ML;
  
    
    u0 = u1;
    u1 = u2;
    u2 = u3;
    u3 = u4; 
   
    uSOR(:,n) = u4([bdaT.EdgeSOR;bdaT.SOR2node]);
    if meshT.solType == 'i'
          uDir(:,n) = u4(bdaT.Dir); 
    end
    if  displayOption.displayF 
        ps = figure(1);
        showsolution(meshT.node,meshT.elem,u4(1:meshT.nC))
        axis(displayOption.axis);
        colorbar
        caxis(displayOption.caxis);
             view(displayOption.viewAngle)
        pause(0.01)
    end
end

if meshT.solType == 'i'
   varargout{1} = uDir;
   varargout{2} = u3;
   varargout{3} = u2;
end

end