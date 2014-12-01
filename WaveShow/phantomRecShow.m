function phantomRecShow
maxRefine = 1;
dx = 2.0e-3;
recStep = 110;
xmin = -0.08;
xmax = 0.08;
ymin = -0.022;
ymax = 0.042;

[node,elem] = squaremesh([xmin,xmax,ymin,ymax], dx); % right triangle uniform mesh
for i = 1 : maxRefine
  [node, elem] = uniformrefine(node,elem);
end


%-----------Split each p2 element into 6 small p1 elements ------%
bary = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
[elem2dof,edge,~] = dofP2(elem); % convert to 6 points P2 elem
edgeMid = (node(edge(:,1),:) + node(edge(:,2),:))/2;
nodeP2= [node;edgeMid;bary]; % convert to masslumping nodes

nD = size(node,1) + size(edge,1);
nE = size(elem,1);
elem2dof = [elem2dof,(nD+1:nD+nE)']; % convert to mass lumping p2 elems

% split the p2 elements into p1 triangles
elemR = [elem2dof(:,1),elem2dof(:,6),elem2dof(:,7);...
         elem2dof(:,1),elem2dof(:,7),elem2dof(:,5);...
         elem2dof(:,5),elem2dof(:,7),elem2dof(:,3);...
         elem2dof(:,7),elem2dof(:,4),elem2dof(:,3);...
         elem2dof(:,7),elem2dof(:,2),elem2dof(:,4);...
         elem2dof(:,6),elem2dof(:,2),elem2dof(:,7)];
     
load(strcat('.\WaveRecData\reconstruction_P2_PE', num2str(recStep),'.mat'), 'f');
showsolution(nodeP2,elemR, 1500./sqrt(1+f));
view([0,-90]);
caxis([1400,1600])
colormap gray

end