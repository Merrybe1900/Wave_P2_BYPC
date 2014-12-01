function markedElem = mark(elem,eta,theta,method,varargin)
% MARK mark element.
%
% markedElem = mark(elem,eta,theta) mark a subset of elements by Dorfler
% marking strategy. It returns an array of indices of marked elements
% markedElem such that sum(eta(markedElem)^2) > theta*sum(eta^2).
%
% markedElem = mark(elem,eta,theta,'max') choose markedElem such that
% eta(markedElem) > theta*max(eta).
%
% markedElem = mark(elem,eta,theta,'COARSEN') choose markedElem such that
% eta(markedElem) < theta*max(eta).
%
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details.

NT = size(elem,1); isMark = false(NT,1);
if ~exist('method','var'), method = 'L2'; end  % default marking is L2 based
switch upper(method)
    case 'MAX'
        isMark(eta>theta*max(eta))=1;
    case 'COARSEN'
        isMark(eta<theta*max(eta))=1;
    case 'L2'
        [sortedEta,idx] = sort(eta.^2,'descend'); 
        x = cumsum(sortedEta);
        isMark(idx(x < theta* x(NT))) = 1;
        isMark(idx(1)) = 1;
    case 'TENSOR_UB'
        isMark(eta<=theta)=1;
    case 'TENSOR_LB'
        bary = varargin{1};
        nodelim = varargin{2};
        isMark(eta>=theta)=1;
        isMark((bary(:,1)<=nodelim(1))|bary(:,2)>=nodelim(2)|bary(:,1)>=nodelim(3)) = 0;
    case 'TENSOR_ULB'
        isMark(eta >= theta(1) & eta <=theta(2))=1;
end
if ~isempty(varargin)&& ~strcmp(method,'TENSOR_LB')
    isMark(varargin{1}>=varargin{2}& varargin{1}<=varargin{3})=1;
end
markedElem = uint32(find(isMark==true));