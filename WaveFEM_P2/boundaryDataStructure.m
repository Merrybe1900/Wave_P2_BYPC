function  bdaT = boundaryDataStructure(meshT,varargin)
%************* boundary data structure function ***************%
% Author: Yun Bai
% Date: 21.08.2014
%**************************************************************%

%----- Source boundary -------%
% bdaT.EdgeSOR: source boundary nodes index (sorted) (M+1)*1.
% bdaT.elSOR: source boundary interval length M*1.
% bdaT.SORid: source edge id.   M*1
%-----------------------------%

% locate the source edge
[~,bdaT.EdgeSOR,bdaT.SORid] = setboundaryWave(meshT,'Neumann', varargin{1});

% horizantal sort
Ic =find(meshT.node(bdaT.EdgeSOR(:,2),1) -meshT.node(bdaT.EdgeSOR(:,1),1)<= 0);
bdaT.EdgeSOR(Ic,:) = [bdaT.EdgeSOR(Ic,2),bdaT.EdgeSOR(Ic,1)]; % switch the index when the order is not consistent

% vertical sort
midNode = (meshT.node(bdaT.EdgeSOR(:,1),1) + meshT.node(bdaT.EdgeSOR(:,2),1))*0.5;
[~,idSOR] = sort(midNode);
bdaT.EdgeSOR = bdaT.EdgeSOR(idSOR,:); % recorder the edge vectors
bdaT.SORid = bdaT.SORid(idSOR,:);     % recorder te edge indices

% compute the edge length
ve = meshT.node(bdaT.EdgeSOR(:,1),:) - meshT.node(bdaT.EdgeSOR(:,2),:);
bdaT.elSOR = sqrt(sum(ve.^2,2));

% store only non-overlap edge nodes
bdaT.EdgeSOR = [bdaT.EdgeSOR(:,1);bdaT.EdgeSOR(end,2)];
bdaT.lenSOR = length(bdaT.EdgeSOR);
bdaT.SOR2node = bdaT.SORid + meshT.nC;

% post check
len = meshT.node(bdaT.EdgeSOR(2:end),1)-meshT.node(bdaT.EdgeSOR(1:end-1));
if max(len)-min(len) > 1e-10
  ('Warning: the emitting boudary is not sorted or lengh is not uniform'); 
end

%----- Absorbing boundary 'ABC'-------%
% bdaT.EdgeABC: source boundary nodes index (sorted). cell 3*1.
% bdaT.elSOR: sourse boundary interval length.  cell 3*1.
% bdaT.ABCid: absorbing edge id. cell 3*1
%-------------------------------------%

nABC = meshT.kABC;
if  meshT.kABC == 3;  % 3 absorbing boundaries
   
    bdaT.EdgeABC = cell(nABC,1);
    bdaT.ABCid = cell(nABC,1);
    bdaT.elABC = cell(nABC,1);
    bdaT.NS = zeros(nABC,1);

    sort_order = {'ascend','descend','descend'};
    sort_id = [2,1,2];

    for k = 1:meshT.kABC

        % locate the kth ABC edge
        [~,bdEdgeT,bdaT.ABCid{k}] = setboundaryWave(meshT,'ABC', varargin{k+1});

        % horizantal sort
        if strcmp(sort_order{k},'ascend')
           Ic =find(meshT.node(bdEdgeT(:,2),sort_id(k)) -meshT.node(bdEdgeT(:,1),sort_id(k))<= 0);
        else
           Ic =find(meshT.node(bdEdgeT(:,2),sort_id(k)) -meshT.node(bdEdgeT(:,1),sort_id(k))>= 0);
        end
        bdEdgeT(Ic,:) = [bdEdgeT(Ic,2),bdEdgeT(Ic,1)];


        midNode = (meshT.node(bdEdgeT(:,1),sort_id(k)) + meshT.node(bdEdgeT(:,2),sort_id(k)))*0.5;
        [~,idt] = sort(midNode,sort_order{k});
        bdaT.ABCid{k} = bdaT.ABCid{k}(idt,:);
        bdEdgeT = bdEdgeT(idt,:);
        bdaT.EdgeABC{k} = [bdEdgeT(:,1);bdEdgeT(end,2)];
        bdaT.NS(k) = length(bdaT.EdgeABC{k});
        bdaT.elABC{k} = sqrt((meshT.node(bdEdgeT(:,1),1)-meshT.node(bdEdgeT(:,2),1)).^2 ...
                  + (meshT.node(bdEdgeT(:,1),2)-meshT.node(bdEdgeT(:,2),2)).^2);

        len = meshT.node(bdaT.EdgeABC{k}(2:end),sort_id(k))-meshT.node(bdaT.EdgeABC{k}(1:end-1),sort_id(k));
        if max(len)-min(len) > 1e-10
           ('Warning: the absorbing boudary is not sorted or lengh is not uniform'); 
        end

    end

    bdaT.bd2node =[bdaT.EdgeABC{1};bdaT.ABCid{1}+meshT.nC; ...
                  bdaT.EdgeABC{2}(2:end);bdaT.ABCid{2}+meshT.nC;...
                  bdaT.EdgeABC{3}(2:end);bdaT.ABCid{3}+meshT.nC];
else
    meshT.kABC = 2;  % 3 absorbing boundaries
    bdaT.EdgeABC = cell(nABC,1);
    bdaT.ABCid = cell(nABC,1);
    bdaT.elABC = cell(nABC,1);
    bdaT.NS = zeros(nABC,1);

    sort_order = {'ascend','descend'};
    sort_id = [2,2];

    for k = 1:2
        % locate the kth ABC edge
        if k == 1
           [~,bdEdgeT,bdaT.ABCid{k}] = setboundaryWave(meshT,'ABC', varargin{k+1});
        elseif k ==2
           [~,bdEdgeT,bdaT.ABCid{k}] = setboundaryWave(meshT,'ABC', varargin{k+2});
        end
        % horizantal sort
        if strcmp(sort_order{k},'ascend')
           Ic =find(meshT.node(bdEdgeT(:,2),sort_id(k)) -meshT.node(bdEdgeT(:,1),sort_id(k))<= 0);
        else
           Ic =find(meshT.node(bdEdgeT(:,2),sort_id(k)) -meshT.node(bdEdgeT(:,1),sort_id(k))>= 0);
        end
        bdEdgeT(Ic,:) = [bdEdgeT(Ic,2),bdEdgeT(Ic,1)];


        midNode = (meshT.node(bdEdgeT(:,1),sort_id(k)) + meshT.node(bdEdgeT(:,2),sort_id(k)))*0.5;
        [~,idt] = sort(midNode,sort_order{k});
        bdaT.ABCid{k} = bdaT.ABCid{k}(idt,:);
        bdEdgeT = bdEdgeT(idt,:);
        bdaT.EdgeABC{k} = [bdEdgeT(:,1);bdEdgeT(end,2)];
        bdaT.NS(k) = length(bdaT.EdgeABC{k});
        bdaT.elABC{k} = sqrt((meshT.node(bdEdgeT(:,1),1)-meshT.node(bdEdgeT(:,2),1)).^2 ...
                  + (meshT.node(bdEdgeT(:,1),2)-meshT.node(bdEdgeT(:,2),2)).^2);

        len = meshT.node(bdaT.EdgeABC{k}(2:end),sort_id(k))-meshT.node(bdaT.EdgeABC{k}(1:end-1),sort_id(k));
        if max(len)-min(len) > 1e-10
           ('Warning: the absorbing boudary is not sorted or lengh is not uniform'); 
        end

    end

    bdaT.bd2node =[bdaT.EdgeABC{1};bdaT.ABCid{1}+meshT.nC; ...
                   bdaT.EdgeABC{2};bdaT.ABCid{2}+meshT.nC];
    
end

% [vt,Ib] = sort(meshT.node(:,2));
% Ic = find(vt<=-0.015);
% notUpdate = Ib(Ic);
% Update = setdiff((1:meshT.nC),notUpdate); 

% ---- Dirichlet nodes ------%
% bdaT.Dir: Dirichlet nodes indices
% bdaT.freenode: non Dirichlet nodes indices 
if meshT.solType == 'i'
    [~,bdaT.Dir, bdaT.DIRid] = setboundaryWave(meshT,'Dirichlet','all');
    bdaT.Dir = [unique(bdaT.Dir);bdaT.DIRid + meshT.nC];
    bdaT.freenode = setdiff(1:meshT.Ndof,bdaT.Dir);
end

end


