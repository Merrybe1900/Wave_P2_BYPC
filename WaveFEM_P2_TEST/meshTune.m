function [mesh,dt,Nt,Ke,Kn] = meshTune(mesh,ref,Nt,dt,soundSpeedWater)

[mesh.elem,~] = fixorientation(mesh.node,mesh.elem);
%% Refinement
if ref.maxRefine > 0
   for i = 1:ref.maxRefine
      [mesh.node, mesh.elem] = uniformrefine(mesh.node,mesh.elem);
      dt = dt*0.5;
      Nt = Nt*2;
   end
end


if ~ref.isUniform
    soundCri = [1560, 1610];
% ----- Obsorbing boundary refinement
%    nE = size(elem,1);
%    ismark = false(nE,1);
%    bary = (node(elem(:,1),:) + node(elem(:,2),:)  + node(elem(:,3),:))/3;
%    ismark((bary(:,1)<=xmin + 2*dxR)|bary(:,2)<= ymin + 2*dxR|bary(:,2)>= ymax - 2*dxR|bary(:,1)>= xmax-2*dxR) = true;
%    bdElem= find(ismark);
%    [node,elem] = bisect(node,elem,bdElem,[]);

    for i = 1:mesh.nRefine
        bary = (mesh.node(mesh.elem(:,1),:) + mesh.node(mesh.elem(:,2),:)  + mesh.node(mesh.elem(:,3),:))/3;
        K = SoundSpeedPhantomBY_ribs(bary,ref.rib_opt);
        if i ==1
          % markedElem1 = mark(mesh.elem,K,soundCriUB,'TENSOR_UB');
           [eta,~] = estimatorAdaptive(mesh.node,mesh.elem,K);
           markedElem2 = mark(mesh.elem,eta,0.015,'MAX');
           markedElem1 = mark(mesh.elem,K,soundCri,'TENSOR_ULB');
           markedElem = union(markedElem1,markedElem2);
        else
        % mesh refinement
            markedElem1 = mark(mesh.elem,K,soundCri,'TENSOR_ULB');
           % markedElem2 = mark(mesh.elem,K,soundCriUB,'TENSOR_UB');
            [eta,~] = estimatorAdaptive(mesh.node,mesh.elem,K);
            markedElem2 = mark(mesh.elem,eta,0.015,'MAX');
            markedElem = union(markedElem1,markedElem2);

        end
        markedElem = checkMarkedElem(bary,dxR,markedElem,1);
        [mesh.node,mesh.elem] = bisect(mesh.node,mesh.elem,markedElem,[]);
        dt = dt*0.5;
        Nt = Nt*2;
    end
end

bary = (mesh.node(mesh.elem(:,1),:) + mesh.node(mesh.elem(:,2),:)  + mesh.node(mesh.elem(:,3),:))/3;
if ref.isHomo % for homogeneous test %
   K0 = soundSpeedWater;
   Kn = 1./(K0.^2)*ones(length(mesh.node),1);
   Ke = 1./(K0.^2)*ones(length(bary),1);
else
   Ke =  SoundSpeedPhantomBY_ribs(bary,ref.rib_opt);
   Ke = 1./(Ke.^2);
   Kn =  SoundSpeedPhantomBY_ribs(mesh.node,ref.rib_opt);
   Kn = 1./(Kn.^2);
end

mesh.nC = size(mesh.node,1);
[mesh.elem,~] = fixorientation(mesh.node,mesh.elem);
sf = strcat('.\MeshData\Reference_mesh_',num2str(mesh.nC),'.mat');
save(sf,'mesh');


end