function WaveInversionMain_P2
%************ Main function of wave inversion simulation ****************%
% Author: Yun Bai
% Date: 21.08.2014
%************************************************************************%

%--- Initial input: User's setting ---%
% all the user modifications are made in this function
[signal,ref,displayOption,mesh,source,Nt,dt,maxCount,...
         expS,expAB1,expAB2,expAB3,sf_Rec,sf_Err,sf_Mesure2,...
                          ini_file,second_ini_file] = RunningExampleInitial;

%--- Received signal---%
% Based on the true phantom to compute measurements on the emitting boundary

if ref.isReceivedSignal
   recievedSignalP2(ref,signal,displayOption,source,Nt,dt,expS,expAB1,expAB2,expAB3)
end
    
%--- Inversion -----%
mesh.solType = 'i';
[mesh.node,mesh.elem,dt,Nt,~,~] = meshRefine(mesh.node,mesh.elem,mesh.solType,ref,Nt,dt);
[mesh.elem2dof,edge,~] = dofP2(mesh.elem); % convert to P2 elem

% compute mesh parameters of reference mesh
mesh.nC = size(mesh.node,1);  
mesh.nE = size(mesh.elem,1); 
nEd = size(edge,1);

if mesh.quadorder == 3.5
   mesh.Ndof = mesh.nC + nEd + mesh.nE; % total DOF
elseif mesh.quadorder == 2
   mesh.Ndof = mesh.nC + nEd;
end

% boundary data structure
bda = boundaryDataStructure(mesh,expS,expAB1,expAB2,expAB3);
IterRec = length(maxCount);  
costF = cell(IterRec,1);
%++++++++ try to introduce 'weight' on line search ++++++++++++%
%The goal of this weight is to have highest weight in the middle and lower
%on the both ends
NodeSup = [mesh.node(:,2);...
          (mesh.node(edge(:,1),2) + mesh.node(edge(:,2),2))*0.5;
          (mesh.node(mesh.elem(:,1),2)...
         + mesh.node(mesh.elem(:,2),2)...
         + mesh.node(mesh.elem(:,3),2))*1/3];
     
mesh.NodeWeight = abs((NodeSup-mesh.ymin)/mesh.rangeY-0.5)*8+1;
mesh.NodeWeight = 1./mesh.NodeWeight*5;
clear NodeSup
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%


for ir = 1:IterRec
    costF{ir} = zeros(maxCount(ir),1);
    load(signal.sf_Mesure,'uSOR');
    
    if ir == 1
       if ~isempty(ini_file)
           load(ini_file,'f'); 
       else
           f = zeros(mesh.Ndof,1);
       end
    else
       load(second_ini_file,'f');
    end
    
    count = 1;
        
        disp('****Inverse reconstruction***')
        if mesh.quadorder == 3.5 %  masslumping solver
            if size(mesh.elem2dof,2) == 6
               mesh.elem2dof = [mesh.elem2dof,(1+nEd+mesh.nC:mesh.nC + nEd + mesh.nE)'];
            end
      
           Kn = mesh.K0^2*(1+f);
           [A,As,area] = stiffMatrixWaveP2(mesh,bda,Kn); % stiffness matrix
           h = 0;

           while count <= maxCount(ir)
                 disp(count);
                 [f, h, costF{ir}(count)] = inverseReconstruction(mesh,bda,source,signal,A,As,area,...
                                          f,uSOR,Nt,dt,displayOption,count,h,sf_Mesure2,ir);

                 if ~mod(count,10)
    %                  showsolution(mesh.node, mesh.elem, 1500./sqrt(1+f(1:mesh.nC)));
    %                  view([0,-90]);
    %                  pause(0.05);
                     countT = count;
                     if ir == 2
                        countT = count + maxCount(1); 
                     end

                     sf = strcat(sf_Rec,num2str(countT),'.mat');
                     save(sf,'f');
                 end
                 count = count + 1;

           end

        elseif mesh.quadorder == 2  %  classical solver

           [A,As,area] = stiffMatrixWaveP2(mesh,bda); %stiffness matrix
           while count <= maxCount
                 [f,h,costF{ir}(count)] = inverseReconstruction(mesh,bda,source,signalt,A,As,area,...
                                         f,uSOR,Nt,dt,displayOption,count,h,sf_Mesure2);

                 count = count + 1;

                 if ~mod(count,10)
                    sf = strcat(sf_Rec,num2str(count),'.mat');
                    save(sf,'f');
                    showsolution(mesh.node, mesh.elem, 1500./sqrt(1+f(1:mesh.nC)));
                    view([0,-90]);
                    pause(0.05);
                 end
           end

        end
end
%save the cost function.
save(sf_Err, 'costF'); 
end