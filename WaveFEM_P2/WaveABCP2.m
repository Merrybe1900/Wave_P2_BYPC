function [phi,phid,phidd,Aux] =  WaveABCP2(u,ut,phi,phid,phidd,Ms,As,NS,bd2node,Kt,dt,...
                                 Ndof,kABC)

% phi is a PN*1 array where N is number of nodes on each boundary
% bdh is the array of edge interval on each boundary
% Assume P = 2 (3rd obsorbing B.C.)
% Assume the corner functions have zeros Neumann condition + du/dt(at corner)

 
%% ------------------ Absorbing Edge i ------------------------------%
parfor i = 1:kABC
    N = NS(i);
    % Form the boundary FEM mass and stiffness matrix (P = 2)
    M = [4*Ms{i},sparse(2*N-1,2*N-1);sparse(2*N-1,2*N-1),4*Ms{i}];
    A = [As{i},As{i};As{i},2*As{i}];
    F = [2*As{i}*u{i};zeros(2*N-1,1)];
    
    F(1) = F(1) + (ut{i}(1))*Kt{i}(1);
    F(N) = F(N) + (ut{i}(2))*Kt{i}(2);

    d = (M/dt/dt + 0.25*A)\((-F - A*(phi{i}+ phid{i}*dt+ 0.25*dt^2*phidd{i}))/dt/dt);
    d1 = phid{i} + 0.5*dt*(d+phidd{i});
    phi{i} = phi{i} + dt*phid{i} + 0.25*dt^2*(phidd{i}+d);
    phidd{i} = d;
    phid{i} = d1;    
end



%% ---------------- Mapping to the main index -----------------%%
Aux = zeros(Ndof,1);

Aux(bd2node) = [phid{1}(1:NS(1)-1); ...
              (phid{1}(NS(1))+ phid{2}(1))*0.5;
               phid{1}((NS(1)+1):(2*NS(1)-1)); ...
               phid{2}(2:NS(2)-1);...
              (phid{2}(NS(2))+ phid{3}(1))*0.5;...
               phid{2}((NS(2)+1):(2*NS(2)-1));
               phid{3}(2:(2*NS(3)-1))];
           
end