function [bdNode,bdEdge,isBdNode,isBdElem] = findboundary(elem,bdFlag,varargin)
%% FINDBOUNDARY finds the boundary of a mesh
%
%   [bdNode,bdEdge,isBdNode] = findboundary(elem) finds boundary nodes and
%   edges for a 2-dimensional mesh. Only the topological structure of
%   the mesh is needed. The boundary edges are consistent with the
%   orientation of triangles, i.e., counterclockwise.
%  
%   [bdNode,bdEdge,isBdNode] = findboundary(elem,bdFlag) finds Dirichlet
%   boundary nodes and Neumann edges.
% 
%   See also findboundary3, setboundary
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = max(elem(:));    
nv = size(elem,2);
if nv == 3 % triangle
    totalEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
elseif nv == 4
    totalEdge = [elem(:,[1,2]); elem(:,[2,3]); elem(:,[3,4]); elem(:,[4 1])];
end    
if exist('bdFlag','var') && ~isempty(bdFlag)
    Dirichlet = totalEdge((bdFlag(:) == 1),:);
    isBdNode = false(N,1); 
    isBdNode(Dirichlet(:)) = true;
    bdNode = find(isBdNode);
    bdEdge = totalEdge((bdFlag(:) == 2) | (bdFlag(:) == 3),:);
else
    totalEdge = sort(totalEdge,2);
    [i,j,s] = find(sparse(double(totalEdge(:,1)),double(totalEdge(:,2)),1));
    bdEdge = [i(s==1),j(s==1)];
    isBdNode = false(N,1); 
    isBdNode(bdEdge) = true;
    bdNode = find(isBdNode);
end
isBdElem = isBdNode(elem(:,1)) | isBdNode(elem(:,2)) | isBdNode(elem(:,3));
end

function dh = minEdgeSize(Edge,node)

dh = (node(Edge(:,1),1) - node(Edge(:,2),1)).^2 ...
   +  (node(Edge(:,1),2) - node(Edge(:,2),2)).^2;
dh= sqrt(min(dh));

end