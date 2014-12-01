function [bdNode,bdFace,isBdNode] = findboundary3(elem,bdFlag)
%% FINDBOUNDARY3 find boundary faces and nodes
%
% [bdNode,bdFace,isBdNode] = findboundary3(elem) find the boundary faces
% and nodes of a 3-dimensional mesh. Note only topological structure of
% the mesh is used. The boundary edges are consistent with the
% orientation of triangles, i.e., counterclockwise.
%
% [bdNode,bdFace,isBdNode] = findboundary3(elem,bdFlag) finds Dirichlet
% boundary nodes and Neumann edges.
% 
% See also findboundary, showboundary3, setboundary3.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = max(elem(:));    
nv = size(elem,2);
if nv == 4
    allFace = [elem(:,[2 4 3]);elem(:,[1 3 4]);elem(:,[1 4 2]);elem(:,[1 2 3])];
    nf = 4;
elseif nv == 8
    allFace = [elem(:,[1 4 3 2]); elem(:,[1 2 6 5]); elem(:,[5 6 7 8]);...
               elem(:,[8 7 3 4]); elem(:,[4 1 5 8]); elem(:,[2 3 7 6])];
   nf = 6;
end
if exist('bdFlag','var')
    Dirichlet = allFace((bdFlag(:) == 1),:);
    isBdNode = false(N,1); 
    isBdNode(Dirichlet(:)) = true;
    bdNode = find(isBdNode);
    bdFace = allFace((bdFlag(:) == 2) | (bdFlag(:) == 3),:);
else
    matlabversion = version;
    if str2double(matlabversion(end-5:end-2)) > 2012
        [face, i2, j] = unique(sort(allFace,2),'rows','legacy');
    else
        [face, i2, j] = unique(sort(allFace,2),'rows');
    end
    NT = size(elem,1);
    i1(j(nf*NT:-1:1)) = nf*NT:-1:1; 
    i1 = i1';
    bdFace = face((i1 == i2),:);
    isBdNode(bdFace) = true;
    bdNode = find(isBdNode);
end