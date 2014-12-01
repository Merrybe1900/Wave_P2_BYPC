function [uSOR, varargout] = waveForwardSolverP2(meshT,bdaT,source,signal,ks,...
                                        A,As,M,Me,Ms,Nt,dt,displayOption,Kn)
%***************** Forward wave classcial function *****************%
% Author: Yun Bai
% Date: 21.08.2014
%********************************************************************%
% Note: For the time being, the auxiliary absorbing boundary function is 
%       not working properly, therefore disabled.
%********************************************************************%

% Initial zero conditions
u0 = zeros(meshT.Ndof,1);
u = zeros(meshT.Ndof,1);

% % Initial absorbing boundary condition
% phi = cell(meshT.kABC,1);
% phid = phi;
% phidd = phi;
% for i = 1:meshT.kABC
%     phi{i} = zeros(2*(2*bdaT.NS(i)-1),1);
%     phid{i} = zeros(2*(2*bdaT.NS(i)-1),1);
%     phidd{i} = zeros(2*(2*bdaT.NS(i)-1),1);
% end
% 
% us = cell(meshT.kABC,1);
% ut = cell(meshT.kABC,1);
% Kt = cell(meshT.kABC,1);

%--- Classcial P2 FEM ---% 
[lambda,weight] = quadpts1(4);
pxy1 = lambda(1,1)*meshT.node(bdaT.EdgeSOR(1:end-1),:) ...
     + lambda(1,2)*meshT.node(bdaT.EdgeSOR(2:end),:);
pxy2 = lambda(2,1)*meshT.node(bdaT.EdgeSOR(1:end-1),:) ...
     + lambda(2,2)*meshT.node(bdaT.EdgeSOR(2:end),:);
pxy3 = lambda(3,1)*meshT.node(bdaT.EdgeSOR(1:end-1),:) ...
     + lambda(3,2)*meshT.node(bdaT.EdgeSOR(2:end),:);
psi = zeros(3,3);
psi(:,1) = (2*lambda(:,1)-1).*lambda(:,1);
psi(:,2) = (2*lambda(:,2)-1).*lambda(:,2);
psi(:,3) = 4*lambda(:,1).*lambda(:,2);
bIn = zeros(bdaT.lenSOR-1,3);

uSOR = zeros(length(bdaT.EdgeSOR)+length(bdaT.SOR2node),Nt);
pSORt = signal.pSOR(ks,:);

if meshT.solType == 'i'
   uDir = zeros(length(bdaT.Dir),Nt); 
end
    
for n=1:Nt  
    
    uIn1 = source(pxy1,pSORt,n*dt,signal);
    uIn2 = source(pxy2,pSORt,n*dt,signal);
    uIn3 = source(pxy3,pSORt,n*dt,signal);
          
  
    bIn(:,1) = (weight(1)*psi(1,1)*uIn1 + weight(2)*psi(2,1)*uIn2...
             + weight(3)*psi(3,1)*uIn3).*bdaT.elSOR;
    bIn(:,2) = (weight(1)*psi(1,2)*uIn1 + weight(2)*psi(2,2)*uIn2...
             + weight(3)*psi(3,2)*uIn3).*bdaT.elSOR;
    bIn(:,3) = (weight(1)*psi(1,3)*uIn1 + weight(2)*psi(2,3)*uIn2...
             + weight(3)*psi(3,3)*uIn3).*bdaT.elSOR;

    b = accumarray([bdaT.EdgeSOR(1:end-1);bdaT.EdgeSOR(2:end);bdaT.SOR2node],bIn(:),[meshT.Ndof 1]); 
     
%     % absorbing boundary condition
%     for i = 1:meshT.kABC
%         us{i} = [u0(bdaT.EdgeABC{i});u0(bdaT.ABCid{i}+ meshT.nC)];
%         ut{i} =[(u(bdaT.EdgeABC{i}(1))- u0(bdaT.EdgeABC{i}(1)))/dt,...
%                 (u(bdaT.EdgeABC{i}(bdaT.NS(i)))- u0(bdaT.EdgeABC{i}(bdaT.NS(i))))/dt];
%         Kt{i} = [sqrt(Kn(bdaT.EdgeABC{i}(1))), sqrt(Kn(bdaT.EdgeABC{i}(bdaT.NS(i))))];
%     end
% 
%     [phi,phid,phidd,Aux] =  WaveABCP2(us,ut,phi,phid,phidd,Ms,As,...
%                                       bdaT.NS,bdaT.bd2node,Kt,...
%                                       dt,meshT.Ndof,meshT.kABC);


    u1 = (M+0.5*Me)\((b-A*u) + M *(2*u-u0) + 0.5*Me*u0);
    u0 = u;
    u = u1;
    uSOR(:,n) = u([bdaT.EdgeSOR;bdaT.SOR2node]);
    
    if meshT.solType == 'i' 
       uDir(:,n) = u(bdaT.Dir);
    end
    
    if  displayOption.displayF 
        ps = figure(1);
        showsolution(meshT.node,meshT.elem,u(1:meshT.nC));
        axis(displayOption.axis);
        colorbar
        caxis(displayOption.caxis);
        view(displayOption.viewAngle)
        pause(0.01)
     end
    
    
end

% for inverse reconstruction, the following data is needed
% 1. solution on all the boundaries at each time step
% 2. solutions over the whole domain for the last two time steps
if meshT.solType == 'i'
   varargout{1} = uDir; 
   varargout{2} = u;
   varargout{3} = u0;
end

end