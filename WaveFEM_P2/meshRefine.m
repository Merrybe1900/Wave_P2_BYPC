function [node,elem,dt,Nt,Ke,Kn] = meshRefine(node,elem,solType,ref,Nt,dt)
%**************** mesh refinement function ***************%
% Author: Yun Bai
% Date: 21.08.2014
%*********************************************************%

%--- make elem indices conter clockwise ---%
[elem,~] = fixorientation(node,elem);

%--- Uniform space refinement ---% 
% time step needs to be refined the same time
if ref.maxRefine > 0
   for i = 1:ref.maxRefine
      [node, elem] = uniformrefine(node,elem);
      dt = dt*0.5;
      Nt = Nt*2;
   end
end

%--- Nonuniform space refinement ---%
if ~ref.isUniform
    soundCri = [1560, 1610];

    for i = 1:mesh.nRefine
        bary = (node(elem(:,1),:) + node(elem(:,2),:)  + node(elem(:,3),:))/3;
        K = SoundSpeedPhantomBY_ribs(bary,ref.rib_opt);
        if i ==1
          % markedElem1 = mark(elem,K,soundCriUB,'TENSOR_UB');
           [eta,~] = estimatorAdaptive(node,elem,K);
           markedElem2 = mark(elem,eta,0.015,'MAX');
           markedElem1 = mark(elem,K,soundCri,'TENSOR_ULB');
           markedElem = union(markedElem1,markedElem2);
        else
            markedElem1 = mark(elem,K,soundCri,'TENSOR_ULB');
           % markedElem2 = mark(elem,K,soundCriUB,'TENSOR_UB');
            [eta,~] = estimatorAdaptive(node,elem,K);
            markedElem2 = mark(elem,eta,0.015,'MAX');
            markedElem = union(markedElem1,markedElem2);
         end
        markedElem = checkMarkedElem(bary,dxR,markedElem,1);
        [node,elem] = bisect(node,elem,markedElem,[]);
        dt = dt*0.5;
        Nt = Nt*2;
    end
    
end

%--- read the reference soundspeed: 1/c^2 ---%
bary = (node(elem(:,1),:) + node(elem(:,2),:)  + node(elem(:,3),:))/3;
if ~ ref.ishomo
    if ref.testType == 6
        Ke =  SoundSpeedPhantomBY_ribs(bary,ref.rib_opt);
        Ke = 1./(Ke.^2);
        Kn =  SoundSpeedPhantomBY_ribs(node,ref.rib_opt);
        Kn = 1./(Kn.^2);
    elseif ref.testType == 7
        Ke =  SoundSpeedPhantomBYsmall_ribs(bary,ref.rib_opt);
        Ke = 1./(Ke.^2);
        Kn =  SoundSpeedPhantomBYsmall_ribs(node,ref.rib_opt);
        Kn = 1./(Kn.^2);
    end
else
%     Ke = 1/(1500^2)*ones(length(bary),1);
%     Kn = 1/(1500^2)*ones(length(node),1);
    Ke = testPhantom(bary,0.16/0.064,ref.testType);
    Ke = 1./(Ke.^2);
    Kn = testPhantom(node,0.16/0.064,ref.testType);
    Kn = 1./(Kn.^2);
end

if solType == 'e'
    sf = strcat(ref.sf_Mesh,'Reference_mesh_',num2str(size(node,1)),'.mat');
    save(sf,'node','elem','-v7.3');
end

end