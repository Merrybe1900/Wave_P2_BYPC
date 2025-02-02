function [x,info,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(A,b,elem,option,varargin)
%% MG multigrid-type solvers
%
%   x = MG(A,b,elem) attempts to solve the system of linear equations A*x =
%   b for x using multigrid type solver. To acheive multigrid efficiency,
%   the mesh must be of bisection type. Inside mg, an coarsening algorithm
%   on bisection grids is applied. See <a href="matlab:ifem coarsendoc">doc coarsen</a> for the coarsening algorithm. 
%
%   The method is designed for the system from several finite element
%   descritzations of elliptic equations on a grid whose topology
%   is given by the array elem. See <a href="matlab:ifem meshdoc">doc mesh</a> for the data structure. 
%
%   x = MG(A,b,elem) works for 2-D piecewise constant or linear finite
%   element by default. For other elements, more mesh structure should be
%   provided in varargin. For 3-D elements, HB information is needed, and
%   should be listed as the first parameter in varargin. Here is a table of
%   possible elements.
%   - mg(A,b,elem)                  2-D linear P1 element
%   - mg(A,b,elem,option,edge)          2-D quadratic P2 element
%   - mg(A,b,elem,option,edge,'CR')     2-D non-conforming CR P1 element
%   - mg(A,b,elem,option,HB)            3-D linear P1 element
%   - mg(A,b,elem,option,HB,edge)       3-D quadratic P2 element
%   - mg(A,b,elem,option,HB,face,'CR')  3-D non-conforming CR P1 element
%
%   x = MG(A,b,elem,options) specifies options in the following list.
%   - option.x0: the initial guess. Default setting x0 = 0.
%   - option.tol: the tolerance of the convergence. Default setting 1e-8.
%   - option.maxIt: the maximum number of iterations. Default setting 200.
%   - option.N0: the size of the coarest grid. Default setting 500.
%   - option.mu-1: smoothing steps
%   - option.coarsegridsolver: solver used in the coarest grid. Default
%     setting: direct solver.
%   - option.solver: various cycles and Krylov space methods
%       * 'NO'     only setup the transfer matrix
%       * 'Vcycle'      V-cycle MultiGrid Method
%       * 'Wcycle'      W-cycle MultiGrid Method
%       * 'Fcycle'      Full cycle Multigrid Method
%       * 'cg'     MG-Preconditioned Conjugate Gradient
%       * 'minres' MG-Preconditioned Minimal Residual Method
%       * 'gmres'  MG-Preconditioned Generalized Minimal Residual Method
%       * 'bicg'   MG-Preconditioned BiConjugate Gradient Method
%       * 'bicgstable' MG-Preconditioned BiConjugate Gradient Stabilized Method
%       * 'bicgstable1' MG-Preconditioned BiConjugate Gradient Stabilized Method
%       The default setting is 'cg' which works well for SPD matrices. For
%       non-symmetric matrices, try 'gmres' and for symmetric but indefinite
%       matrices, try 'minres' or 'bicg' sequences.
%       The string option.solver is not case sensitive.
%   - option.preconditioner:  multilevel preconditioners including:
%       * 'V'   V-cycle MultiGrid used as a Preconditioner
%       * 'W'   W-cycle MultiGrid used as a Preconditioner
%       * 'F'   Full cycle Multigrid used as a Preconditioner
%       * 'bpx' BPX-Preconditioner
%   - option.printlevel: the level of screen print out
%       * 0: no output
%       * 1: name of solver and convergence information (step, err, time)
%       * 2: convergence history (err in each iteration step)
%
%   [x,info] = MG(A,b,elem) also returns information of the solver
%   - info.flag:
%       * 0: mg converged to the desired tolerance tol within maxIt iterations
%       * 1: mg iterated maxIt times but did not converge.
%       * 2: direct solver
%   - info.itStep: the iteration number at which x was computed.
%   - info.time: the cpu time to get x
%   - info.err: the approximate relative error in the energy norm in
%   err(:,1) and the relative residual norm(b-A*x)/norm(b) in err(:,2). If
%   flag is 0, then max(err(end,:)) <= tol.
%   - info.stopErr: the error when iteration stops
%
%   Example:
%
%   [node,elem] = squaremesh([0 1 0 1],0.5);
%   for k = 1:8
%     [node,elem] = uniformrefine(node,elem);
%   end
%   pde.f = inline('p(:,1).*p(:,2)','p');
%   pde.g_D = inline('zeros(size(p,1),1)','p');
%   option.solver = 'notsolve';
%   [u,Du,eqn] = Poisson(node,elem,pde,[],option);
%   fprintf('\n Number of unknowns: %8.0u\n',length(eqn.b))
%   tic; display('Direct solver'); u = eqn.A\eqn.b; toc;
%   tic; x = mg(eqn.A,eqn.b,elem); toc;
%   format shorte
%   fprintf('Difference between direct and mg solvers %0.2g \n',norm(u-x));
%
% See also mgMaxwell, dmg
%
% Documentation in Help browser <a href="matlab:ifem mgdoc">ifem mgdoc</a>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

tic;

%% Size of systems
Ndof = length(b);                  % number of dof
N = max(elem(:));                  % number of nodes
NT = size(elem,1);                 % number of elements
dim = size(elem,2)-1;
if N > NT       % 2-D quad mesh
    dim = 2;
end

%% Options
% Assign default values to unspecified parameters
if ~exist('option','var'), 
    option = []; 
end
option = mgoptions(option,Ndof);    % parameters
x0 = option.x0; 
N0 = option.N0; 
tol = option.tol;
maxIt = option.solvermaxit; 
mu = option.smoothingstep; 
solver = option.solver; 
preconditioner = option.preconditioner;
coarsegridsolver = option.coarsegridsolver; 
printlevel = option.printlevel; 
setupflag = option.setupflag;

if setupflag == true
    isFixDof = [];
    if isfield(option,'freeDof') % freeDof is given
       isFreeDof = false(N,1);   
       isFreeDof(option.freeDof) = true;
       NA = size(A,1);
       if NA > length(option.freeDof)
           isFixDof = true(NA,1);
           isFixDof(isFreeDof) = false;
       end
    else % Find free dof and eliminate isolated dof
        deg = sum(spones(A));  % degree 
        isFreeDof = false(Ndof,1);
        isFreeDof(deg>1) = true;
        isFixDof = ~isFreeDof;
        isFixDof(deg == 0) = false;
    end
    if any(isFixDof) % a bigger matrix is given
        xD = zeros(Ndof,1);
        xD(isFixDof) = b(isFixDof)./diag(A(isFixDof,isFixDof));
        A = A(isFreeDof,isFreeDof);
    end
    if length(b) > sum(isFreeDof)
        b = b(isFreeDof);
    end
    if length(x0) > sum(isFreeDof)
        x0 = x0(isFreeDof);
        option.x0 = x0;
    end
isFreeNode = isFreeDof(1:N); % default choice for P1

%% Additional Hierarchical Structure for Different Elements
if Ndof > N % other than P1 element
    NE = N + NT -1;   % estimate of NE by Euler formula
    if dim == 3
        NF = 2*NT; % estimate of NF
    end
    % get additional data structure 
    if Ndof ~= NT     % not piecewise constant element
        if dim == 2
            if nargin > 4
                edge = varargin{1};
                NE = size(edge,1);
            else
                [tempvar,edge] = dofedge(elem);
            end
        elseif dim == 3
            if nargin > 5 % additional data structure
                if size(varargin{1},2) == 2 % edge 
                    edge = varargin{1};
                    NE = size(edge,1);
                elseif size(varargin{1},2) == 3 % face
                    face = varargin{1};
                    NF = size(face,1);
                end
            else % if no additional edge/face, generate it
                 [tempvar,edge] = dof3edge(elem); %#ok<*ASGLU>
                 [tempvar,face] = dof3face(elem);
                 NE = size(edge,1);
                 NF = size(face,1);
            end
        end
    end
    % transfer operator from P1 to the current element
    isFreeNode = true(N,1); 
    if Ndof == NT % piecewise constant element        
       if dim == 2
%           [bdNode,bdEdge,isBdNode] = findboundary(elem); 
%           isFreeNode = ~isBdNode;
          P1toP0 = sparse([1:NT;1:NT;1:NT]',elem,ones(3*NT,1)/3,NT,N); 
       elseif dim == 3
%           [bdNode,bdFace,isBdNode] = findboundary3(elem); 
%           isFreeNode = ~isBdNode;
          P1toP0 = sparse([1:NT;1:NT;1:NT;1:NT]',elem,ones(4*NT,1)/4,NT,N); 
       end
       auxPro = P1toP0(isFreeDof,1:N);            
       isFreeNode = [];
    end
    if Ndof == N + NE % quadratic element 
        isFreeNode = isFreeDof(1:N); 
        P1toP2 = sparse([(1:N)'; N+(1:NE)'; N+(1:NE)'], ...
                        [(1:N)'; double(edge(:))],...
                        [ones(N,1); 0.5*ones(2*NE,1)]',Ndof,N);
        auxPro = P1toP2(isFreeDof,isFreeNode);        
    end
    if (dim == 2) && (Ndof == N + 2*NE +NT) % cubic element in 2-D
        isFreeNode = isFreeDof(1:N); 
        P1toP3 = sparse([(1:N)'; N+2*(1:NE)'-1; N+2*(1:NE)'-1; N+2*(1:NE)';...
                          N+2*(1:NE)';N+2*NE+(1:NT)';N+2*NE+(1:NT)';N+2*NE+(1:NT)'], ...
                        [(1:N)'; double(edge(:));double(edge(:));elem(:);       ],...
                        [ones(N,1); 2/3*ones(NE,1);1/3*ones(NE,1);...
                         1/3*ones(NE,1);2/3*ones(NE,1);1/3*ones(3*NT,1)]',Ndof,N);
         auxPro = P1toP3(isFreeDof,isFreeNode);
    end
    if (dim == 2) && (Ndof == NE)  % 2-D CR nonconforming element
        isFreeNode(edge(isFixDof,:)) = false;
        P1toCR = sparse([1:NE;1:NE]',double(edge(:)),0.5*ones(2*NE,1),NE,N);
        auxPro = P1toCR(isFreeDof,isFreeNode);        
    end
    if (dim == 3) && (Ndof == NF) % 3-D CR nonconforming element
        isFreeNode(face(isFixDof,:)) = false;
        P1toCR = sparse([1:NF;1:NF;1:NF]',double(face(:)),ones(3*NF,1)/3,NF,N);        
        auxPro = P1toCR(isFreeDof,isFreeNode);        
    end
    if (dim == 2) && (Ndof == NT + NE) % 2-D weak Galerkin element (P0,P0,RT0) element
%         if (nargin > 4) && strcmp(varargin{end},'WG') % 
        fixEdgeDof = find(isFixDof) - NT;
        isFreeNode(edge(fixEdgeDof,:)) = false;
        P1toWG = sparse([repmat((1:NT)',3,1); repmat((NT+1:Ndof)',2,1)], ... % i
                        [elem(:);             double(edge(:))], ...   % j 
                        [ones(3*NT,1)/3;      ones(2*NE,1)/2], NT+NE, N);        
        auxPro = P1toWG(isFreeDof,isFreeNode);        
    end
    if (dim == 3) && (Ndof == NT + NF) % 3-D weak Galerkin element (P0,P0,RT0) element
%         if (nargin > 4) && strcmp(varargin{end},'WG') % 
        fixFaceDof = find(isFixDof) - NT;
        isFreeNode(face(fixFaceDof,:)) = false;
        P1toWG = sparse([repmat((1:NT)',4,1); repmat((NT+1:Ndof)',3,1)], ... % i
                        [elem(:);           double(face(:))], ...   % j 
                        [ones(4*NT,1)/4;    ones(3*NF,1)/3], NT+NF, N);
        auxPro = P1toWG(isFreeDof,isFreeNode);        
    end
    if (nargin > 4) && (dim == 2) && ischar(varargin{1}) && strcmp(varargin{1},'HB') % HB basis element
        isFreeNode = isFreeDof(1:N); 
        P1toPk = speye(Ndof,N);
        auxPro = P1toPk(isFreeDof,isFreeNode);        
    end
end

%% Hierarchical Structure of Mesh
%  hierarchical structure for linear element.
if dim == 2  % 2-D
    [HB, NL, level] = HBstructure(elem,N0); 
end
if dim == 3  % 3-D
    if nargin > 4
        HBmesh = varargin{end};
    end
    [HB, NL, level] = HBstructure3(elem,HBmesh,N0);
end

%% No coarsening or coarsened nodes is small
nnzlevel = find(NL);
if Ndof <= N    % linear element
    if (level == 1) || (NL(nnzlevel(end)) - NL(nnzlevel(1)))/N < 0.25
        if Ndof < N0 % small size
            solver = 'DIRECT';
        else
            [x,info] = amg(A,b,option);
            if any(isFixDof)
                xD(isFreeDof) = x;
                x = xD;
            end
            Ai{1} = A; Bi{1} = tril(A); BBi{1} = Bi'; Res = []; Pro = [];
            return
        end
    end
end    

%% Direct solver
if strcmp(solver,'DIRECT')        
    x = A\b;                       % use direct solver                          
    xD(isFreeDof) = x;
    x = xD;
    flag = 2; itStep = 0; err = norm(b-A*x)/norm(b); time = toc;
    if printlevel >= 1
        fprintf('Direct solver \n')
        fprintf('#dof: %8.0u,  #nnz: %8.0u,  iter: %2.0u,  err = %8.4e,  time = %4.2g s\n',...
                 Ndof, nnz(A), itStep, max(err(end,:)), time)
        if level == 1
            fprintf('No coarsening of the grid. \n')
        end
    end
    info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));    
    Ai{1} = A; Bi{1} = tril(A); BBi{1} = Bi'; Res = []; Pro = [];
    return
end


%% Transfer operators between multilevel meshes for P1 element
% standard prolongation and restriction operator for P1 element
[Pro,Res] = transferoperator(HB,NL,isFreeNode); 
if Ndof > N
    if ~exist('auxPro','var')
        disp('The current element is not supported by mg');
    else
        Pro{level} = auxPro;
        Res{level+1} = Pro{level}'; 
        level = level + 1; % add one more level from P1 to current element
    end
end
clear HB auxPro

%% Matrices in each level
Ai = cell(level,1);
Ai{level} = A;    
for j = level:-1:2
    Ai{j-1} = Res{j}*Ai{j}*Pro{j-1};           % Ac = Res*Af*Pro
    switch option.smoother
        case 'GS'
            Bi{j} = tril(Ai{j});        % Forward Gauss-Seidel   B = D+L
            BBi{j} = triu(Ai{j});       % Backward Gauss-Seidel BB = D+U    
        case 'JAC'
            Bi{j} = spdiags(diag(Ai{j}),0,size(Ai{j},1),size(Ai{j},1));        %
            BBi{j} = Bi{j};       % Jacobi iteration    
    end
    if option.smoothingparameter~=1
        Bi{j} = Bi{j}/option.smoothingparameter;
        BBi{j} = BBi{j}/option.smoothingparameter;
    end
end
% if condest(Ai{1}) > 1e6         % nearly singular
%     Ai{1} = Ai{1} + 1e-3*spdiags(diag(Ai{1}),0,size(Ai{1},1),size(Ai{1},1));
% end
end % end for setup

%% Only set up the transfer operators
if strcmp(solver,'NO') 
    x = x0; flag = 0; itStep = 0; err = 0; time = toc;
    info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));
    return
end

%% No need of set up
if setupflag == false
    Ai = varargin{1};
    Bi = varargin{2};
   BBi = varargin{3};
   Res = varargin{4};
   Pro = varargin{5};
   if isequal(length(varargin),6) 
      isFreeDof = varargin{6};
   else
      isFreeDof = true(Ndof,1);
   end
   level = length(Ai);
   isFixDof = ~isFreeDof;
   xD = zeros(Ndof,1);
   if any(isFixDof)
        xD(isFixDof) = b(isFixDof)./diag(A(isFixDof,isFixDof));
        b = b(isFreeDof);
        A = A(isFreeDof,isFreeDof);
        x0 = x0(isFreeDof);
   end
end

%% Krylov iterative methods use Multigrid-type Preconditioners
% set up preconditioner
if ~strcmp(solver(2:end),'CYCLE')
    switch preconditioner
        case 'V'
            prefunc = @vcycle;
            if printlevel >= 1
                fprintf('Multigrid V-cycle Preconditioner with ')
            end
        case 'W'
            prefunc = @wcycle;
            if printlevel >= 1
                fprintf('Multigrid W-cycle Preconditioner with ')
            end
        case 'F'
            prefunc = @fcycle;
            if printlevel >= 1
                fprintf('Multigrid Full Cycle Preconditioner with ')
            end
        case 'BPX'
            % modify smoother
            for j = level:-1:2
                Di{j} = diag(Ai{j});
            end
            prefunc = @bpx;
            if printlevel >= 1
                fprintf('BPX Preconditioner with ')
            end
    end
end
% initial set up
k = 1; 
x = x0;
r = b - A*x;
nb = norm(b);
err = zeros(maxIt,2);
if nb > eps  % nb is non-zero
    err(1,:) = norm(r)/nb; 
else
    err(1,:) = norm(r);
end
% solvers
switch solver
    case 'CG'
        if printlevel >= 1
            fprintf('Conjugate Gradient Method\n')
        end
        while (max(err(k,:)) > tol) && (k <= maxIt)    
            % compute Br by MG
            Br = prefunc(r);
            % update tau, beta, and p
            rho = Br'*r;  % e'*ABA*e approximates e'*A*e
            if k == 1
                p = Br;
            else
                beta = rho/rho_old;
                p = Br + beta*p;
            end
            % update alpha, x, and r
            Ap = A*p;
            alpha = rho/(Ap'*p);
            r = r - alpha*Ap;
            x = x + alpha*p;
            rho_old = rho;
            k = k + 1;
            % compute err for the stopping criterion
        %     err(k,1) = alpha*sqrt(p'*Ap/(x'*A*x)); % increamental error in energy norm
            err(k,1) = sqrt(abs(rho/(x'*b))); % approximate relative error in energy norm
            err(k,2) = norm(r)/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u,  #nnz: %8.0u, MGCG iter: %2.0u, err = %8.4e\n',...
                         Ndof, nnz(A), k-1, max(err(k,:)));
            end
        end
        err = err(1:k,:);
        itStep = k-1;
    case 'VCYCLE'  
        if printlevel >= 1
            fprintf('Multigrid Vcycle Iteration \n')
        end
        while (max(err(k,:)) > tol) && (k <= maxIt)
            k = k + 1;
            % Step 2: Compute Br by one Vcylce MG
            Br = vcycle(r);
            % Step 3: Correct the solution
            x = x + Br;
            err(k,1) = sqrt(abs(Br'*r/(x'*b))); % approximate relative error in energy norm
            % Step 1: Form residual r
            r = r - A*Br;
            err(k,2) = norm(r)/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u,  #nnz: %8.0u, MG Vcycle iter: %2.0u, err = %8.4e\n',...
                         Ndof, nnz(A), k-1, max(err(k,:)));
            end            
        end
        err = err(1:k,:);
        itStep = k-1;
    case 'WCYCLE'  
        if printlevel >= 1
            fprintf('Multigrid Wcycle Iteration \n')
        end
        while (max(err(k,:)) > tol) && (k <= maxIt)
            k = k + 1;
            % Step 2: Compute Br by one Vcylce MG
            Br = wcycle(r);
            % Step 3: Correct the solution
            x = x + Br;
            err(k,1) = sqrt(abs(Br'*r/(x'*b))); % approximate relative error in energy norm
            % Step 1: Form residual r
            r = r - A*Br;
            err(k,2) = norm(r)/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u,  #nnz: %8.0u, MG Wcycle iter: %2.0u, err = %8.4e\n',...
                         Ndof, nnz(A), k-1, max(err(k,:)));
            end            
        end
        err = err(1:k,:);
        itStep = k-1;
    case 'FCYCLE'  
        if printlevel >= 1
            fprintf('Multigrid Full Cycle Iteration \n')
        end
        while (max(err(k,:)) > tol) && (k <= maxIt)
            k = k + 1;
            % Step 2: Compute Br by one Full cycle MG
            Br = fcycle(r);
            % Step 3: Correct the solution
            x = x + Br;
            err(k,1) = sqrt(abs(Br'*r/(x'*b))); % approximate relative error in energy norm
            % Step 1: Form residual r
            r = r - A*Br;
            err(k,2) = norm(r)/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u, #nnz: %8.0u, MG Fcycle iter: %2.0u, err = %8.4e\n',...
                         Ndof, nnz(A), k-1, max(err(k,:)));
            end            
        end
        err = err(1:k,:);
        itStep = k-1;
    case 'MINRES'
        if printlevel >= 1
            fprintf('Minimum Residual Method \n')
        end
        [x,flag,err,itStep] = minres(A,b,tol,maxIt,prefunc,[],x0);  
    case 'GMRES'
        if printlevel >= 1
            fprintf('General Minimum Residual Method\n')
        end
        if isfield(option,'restart')
            restart = option.restart;
        else
            restart = min(N,10);
        end
        [x,flag,err,itStep] = gmres(A,b,restart,tol,maxIt,prefunc,[],x0);
        itStep = (itStep(1)-1)*restart + itStep(2);
    case 'BICG'
        if printlevel >= 1
            fprintf('BiConjugate Gradient Method\n')
        end
        [x,flag,err,itStep] = bicg(A,b,tol,maxIt,prefunc,[],x0);  
    case 'BICGSTAB'
        if printlevel >= 1
            fprintf('Stablilized BiConjugate Gradient Method\n')
        end
        [x,flag,err,itStep] = bicgstab(A,b,tol,maxIt,prefunc,[],x0);  
    case 'BICGSTAB1'
        if printlevel >= 1
            fprintf('Stablilized BiConjugate Gradient Method\n')
        end
        [x,flag,err,itStep] = bicgstab1(A,b,tol,maxIt,prefunc,[],x0);  
end

%% Modify x to include fix dof
if any(isFixDof)
    xD(isFreeDof) = x;
    x = xD;
end

%% Output
if k > maxIt
    flag = 1;
else
    flag = 0;
end
time = toc;
if printlevel >= 2
    fprintf('#dof: %8.0u, level: %2.0u,   coarse grid %2.0u, #nnz: %8.0u\n',...
              Ndof, level, size(Ai{1},1), nnz(Ai{1}))
end
if printlevel >= 1
    fprintf('#dof: %8.0u,  #nnz: %8.0u, iter: %2.0u,   err = %8.4e,   time = %4.2g s\n',...
                 Ndof, nnz(A), itStep, max(err(end,:)), time)
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions vcycle, wcycle, fcycle, bpx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vcycle MG
    function Br = vcycle(r,J)        % solve equations Ae = r in each level  
    if nargin<=1
        J = level;
    end
    ri = cell(J,1);            % record residual in each level
    ei = cell(J,1);            % record err in each level
    ri{J} = r;
    for i = J:-1:2
        ei{i} = Bi{i}\ri{i};   % pre-smoothing
        for s = 1:mu-1           % extra mu-1 steps smoothing
            ei{i} = ei{i} + Bi{i}\(ri{i}-Ai{i}*ei{i}); 
        end
        ri{i-1} = Res{i}*(ri{i} - Ai{i}*ei{i});
    end
    if strcmp(coarsegridsolver,'direct')
        ei{1} = Ai{1}\ri{1}; % direct solver in the coarest level
    else                         % iterative solver in the coarest level
        D = spdiags(diag(Ai{1}),0,size(Ai{1},1),size(Ai{1},1));
        [ei{1},flag] = pcg(Ai{1},ri{1},1/size(Ai{1},1),1000,D);
    end
    for i = 2:J
        ei{i} = ei{i} + Pro{i-1}*ei{i-1};
        ei{i} = ei{i} + BBi{i}\(ri{i}-Ai{i}*ei{i});
        for s = 1:mu-1
            ei{i} = ei{i} + BBi{i}\(ri{i}-Ai{i}*ei{i}); % post-smoothing
        end
    end
    Br = ei{J};
    end

%% Wcycle MG
    function e = wcycle(r,J)        % solve equations Ae = r in each level  
    if nargin<=1
        J = level;
    end
    if J == 1
        e = Ai{J}\r;   % exact solver in the coaresest grid
        return
    end
    % fine grid pre-smoothing
    e = Bi{J}\r;        % pre-smoothing
    for s = 1:mu-1        % extra mu-1 steps smoothing
        e = e + Bi{J}\(r-Ai{J}*e); 
    end
    % coarse grid correction twice
    rc = Res{J}*(r - Ai{J}*e);
    ec = wcycle(rc,J-1);
    ec = ec + wcycle(rc - Ai{J-1}*ec,J-1);
    e = e + Pro{J-1}*ec;
    % fine grid post-smoothing
    e = e + BBi{J}\(r-Ai{J}*e);
    for s = 1:mu-1
        e = e + BBi{J}\(r-Ai{J}*e); % post-smoothing
    end
    end

%% Fcycle MG
    function Br = fcycle(r)
    ri = cell(level,1);            % record residual in each level
    ei = cell(level,1);            % record err in each level
    ri{level} = r;
    for i = level:-1:2
        ei{i} = vcycle(ri{i},i);   % pre-smoothing
        for s = 1:mu-1               % extra smoothing steps
            ei{i} = ei{i} + vcycle(ri{i}-Ai{i}*ei{i},i); % pre-smoothing
        end
        ri{i-1} = Res{i}*(ri{i} - Ai{i}*ei{i});
    end
    if strcmp(coarsegridsolver,'direct')
        ei{1} = Ai{1}\ri{1};        % direct solver in the coarest level
    else                            % iterative solver in the coarest level
        D = spdiags(diag(Ai{1}),0,size(Ai{1},1),size(Ai{1},1));
        [ei{1},flag] = pcg(Ai{1},ri{1},1/size(Ai{1},1),1000,D);
    end
    for i = 2:level
        ei{i} = ei{i} + Pro{i-1}*ei{i-1};
        ei{i} = ei{i} + vcycle(ri{i}-Ai{i}*ei{i},i);
        for s = 1:mu-1   % extral smoothing
            ei{i} = ei{i} + vcycle(ri{i}-Ai{i}*ei{i},i); % post-smoothing
        end
    end
    Br = ei{level};    
    end

%% BPX preconditioner
    function Br = bpx(r)
    ri = cell(level,1);            % record residual in each level
    ei = cell(level,1);            % record err in each level
    % compute Br by BPX
    ri{level} = r;
    for i = level:-1:2
        ei{i} = ri{i}./Di{i};       % Jacobi smoothing
        ri{i-1} = Res{i}*ri{i};     % restriction of the residual
    end
    if strcmp(coarsegridsolver,'direct')
        ei{1} = Ai{1}\ri{1};         % direct solver in the coarest level
    else                           % iterative solver in the coarest level
        D = spdiags(diag(Ai{1}),0,size(Ai{1},1),size(Ai{1},1));
        [ei{1},flag] = pcg(Ai{1},ri{1},1/size(Ai{1},1),1000,D);
    end
    for i=2:level
        ei{i} = ei{i} + Pro{i-1}*ei{i-1};  % prolongation of the correction
    end
    Br = ei{level};
    end

end
