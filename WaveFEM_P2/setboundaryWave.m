function [bdFlag, varargout] = setboundaryWave(mesh,varargin)
%% SETBOUNDARY set type of boundary edges.
%
%  bdFlag = SETBOUNDARY(node,elem,'Dirichlet') set all boundary edges to
%  Dirichlet type. 
%
%  bdFlag = SETBOUNDARY(node,elem,'Neumann') set all boundary edges to
%  Neumann type. 
%
%  bdFlag = SETBOUNDARY(node,elem,'Robin') set all boundary edges to
%  Robin type (can be replaced by 'ABC'). 
%
%  bdFlag = SETBOUNDARY(node,elem,'Dirichlet','(x==1) | (x==-1)') set
%  Dirichlet boundary condition on x=1 and x=-1. Other edges are
%  homongenous Neumann boundary condition.
%
%  bdFlag = SETBOUNDARY(node,elem,'Dirichlet','(x==1) | ...
%  (x==-1)','Neumann','(y==1) | (y==-1)') set
%  Dirichlet boundary condition on x=1 or x=-1 and Neumann boundary
%  condition on y=1 or y=-1.
%
%  bdFlag = SETBOUNDARY(node,elem,'Dirichlet','(x==1) | ...
%  (x==-1)','Neumann','y==1', 'Robin',' y==-1') set
%  Dirichlet boundary condition on x=1 or x=-1 and Neumann boundary
%  condition on y=1, and Robin boundary condition on y=-1.
%
%  bdFlag = SETBOUNDARY(node,elem,'Dirichlet','all','Neumann','y==1') set
%  Neumann boundary condition on y=1 and others are Dirichlet boundary condition.
%
% Example
%   
%      node = [0,0; 1,0; 1,1; 0,1];
%      elem = [2,3,1; 4,1,3];
%      bdFlag = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');
%      [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
%      [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
%      showmesh(node,elem);
%      allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
%      Dirichlet = allEdge((bdFlag(:) == 1),:);
%      Neumann = allEdge((bdFlag(:) == 2) | (bdFlag(:) == 3),:);
%      findedge(node,Dirichlet,[],'noindex','LineWidth',4,'Color','r');
%      findedge(node,Neumann,[],'noindex','LineWidth',4,'Color','b');
%
% See also setboundary3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Find boundary edges
elem = mesh.elem;
node = mesh.node;
nv = size(elem,2);
if nv == 3 % triangles 
    allEdge = uint32(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
elseif nv == 4 % quadrilateral
    allEdge = uint32(sort([elem(:,[1,2]); elem(:,[2,3]); elem(:,[3,4]); elem(:,[4 1])],2));    
end
ne = nv; % number of edges in one element
Neall = length(allEdge);
matlabversion = version;
if str2double(matlabversion(end-5:end-2)) > 2012
    [edge, i2, j] = unique(allEdge,'rows','legacy'); 
else
    [edge, i2, j] = unique(allEdge,'rows'); 
end
NT = size(elem,1);
i1(j(Neall:-1:1)) = Neall:-1:1; 
i1 = i1';
bdFlag = zeros(Neall,1,'uint8');
bdEdgeidx = i1(i1==i2);

%% Set up boundary edges
nVarargin = size(varargin,2);
if (nVarargin==1)
    bdType = findbdtype(varargin{1});
    bdFlag(bdEdgeidx) = bdType;
end
if (nVarargin>=2)
    for i=1:nVarargin/2
        bdType = findbdtype(varargin{2*i-1});
        expr = varargin{2*i};
        if strcmp(expr,'all')
            bdFlag(bdEdgeidx) = bdType;
        else
           x = (node(allEdge(bdEdgeidx,1),1) + node(allEdge(bdEdgeidx,2),1))/2; %#ok<NASGU>
           y = (node(allEdge(bdEdgeidx,1),2) + node(allEdge(bdEdgeidx,2),2))/2; %#ok<NASGU>
           idx = eval(expr);
           bdFlag(bdEdgeidx(idx)) = bdType;
        end
    end
end
bdFlag = reshape(bdFlag,NT,ne);

if ~strcmp(expr,'all')
    bdEdge = allEdge(bdEdgeidx(idx),:);
    [~,bdEdgeId,~] = intersect(edge,bdEdge,'rows');
else
    bdEdge = allEdge(bdEdgeidx,:);
    [~,bdEdgeId,~] = intersect(edge,bdEdge,'rows');
end

varargout{1} = bdEdge;  % boundary edge nodes indices
varargout{2} = bdEdgeId;   % boundary to edge index to order to locate the edge middle points  


end
%%
function bdType = findbdtype(bdstr)
    switch bdstr
        case 'Dirichlet'
            bdType = 1;
        case 'Neumann'
            bdType = 2;
        case 'Robin'
            bdType = 3;
        case 'ABC' % absorbing boundary condition for wave-type equations
            bdType = 3;
    end
end
