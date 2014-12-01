function recievedSignalP2(ref,signal,displayOption,source,Nt,dt,...
                          expS,expAB1,expAB2,expAB3)
%**************** received signal simulation function ***************%
% Author: Yun Bai
% Date: 21.08.2014
%*********************************************************************%

disp('*** this function is to compute the received signal for measurements *** ')
%--------------------- Input parameter -----------------------------%
% ref: structure contains reference forward setting
% signal: structure, emitting signal parameter
% source: emitting signal
% Nt: total step number; dt: time step size
% expS: string, expression of emitting boudary,e.g. 'y = -0.022'
% expAB1,expAB2,expAB3: three absorbing boudary expressions
%-------------------------------------------------------------------%

tic
ref.solType = 'e';
% refine the initial mesh
[ref.node,ref.elem,dt,Nt,Ke,Kn] = meshRefine(ref.node,ref.elem,...
                                            ref.solType,ref,Nt,dt);
[ref.elem2dof,edge,~] = dofP2(ref.elem); % convert to P2 elem

% compute mesh parameters of reference mesh
ref.nC = size(ref.node,1);  
ref.nE = size(ref.elem,1); 
nEd = size(edge,1);

showsolution(ref.node,ref.elem,1./sqrt(Kn(1:ref.nC)))
view(0,-90)
colorbar
% boundary data structure
bda = boundaryDataStructure(ref,expS,expAB1,expAB2,expAB3);

% Domain stiffness matrix and boundary stiffness matrix % 
if ref.quadorder == 3.5 % masslumping solver
   
   ref.elem2dof = [ref.elem2dof,(1+nEd+ref.nC:ref.nC + nEd + ref.nE)'];
   ref.Ndof = ref.nC + nEd + ref.nE;
   [A,As,area] = stiffMatrixWaveP2(ref,bda,Kn); % stiffness matrix

elseif ref.quadorder == 2 % classical solver
   
   ref.Ndof = ref.nC + nEd;
   [A,As,area] = stiffMatrixWaveP2(ref,bda); %stiffness matrix
   
end

disp(strcat('meshsize is ',num2str(ref.nC), ' DOF = ', num2str(ref.Ndof)));

% Domain mass matrix and boundary mass matrix %
[M,Me,Ms] = massMatrixWaveP2(ref,bda,area,Ke,Kn,dt);

% Forward solver %
uSOR = cell(signal.numS,1);
if ref.quadorder == 3.5 % masslumping 
    parfor ks = 1 : signal.numS
           uSOR{ks} = waveForwardSolverP2massLumping2(ref,bda,source,signal,...
                                    ks,A,As,M,Me,Ms,Nt,dt,displayOption,Kn);
    end
elseif ref.quadorder == 2 % classical
    parfor ks = 1 : signal.numS
           uSOR{ks} = waveForwardSolverP2(ref,bda,source,signal,...
                                    ks,A,As,M,Me,Ms,Nt,dt,displayOption,Kn);
    end
end

tR = toc;
clear bda
disp(strcat('Forward computation takes ', num2str(tR),' seconds'));
save(signal.sf_Mesure,'uSOR','tR','-v7.3');

end