function WaveForwardTest_P2
global xmin xmax ymax

soundSpeedWater = 1500; % [m/s]
K0 = 1/soundSpeedWater; % constant 1/c0
pde = SignalSetup;

%--- signal description ---%
signal.signalDuration = 150e-6; % [s]
signal.samplingFrequency = 5e6; % [Hz]
signal.fc = 150e3; % [Hz]
signal.bw = 150e3;       % [Hz]
signal.sigtype = 'exp';      % 'exp' or 'sinc'
signal.t0 = 2e-6;        % [s]
signal.am = 120;
%pde = SignalSetup;

%-- Computational domain description ---%
xmin = -0.08;
xmax = 0.08;
ymin = -0.022;
ymax = 0.042;
signal.dx = 2.0000e-03;         % initial coarse meshsize
dx_sensor = 0.5*signal.dx;


%-- Reference or Test ---%
errComp = 1;
ref.maxRefine = 1;
ref.isUniform = 1;
ref.rib_opt = 1;
ref.isHomo = 0;  % if ture, homogenous test
% Nref = 2;  % if refinement is nonuniform, then need to indicate the times
             % of nonuniform refinement
             

%-- Display option ---%

displayF = 1;         % Forward solution display
maxZ = 0.08;               % Display range
minZ = -0.03;
displayOption.viewAngle = [0,-90];
displayOption.axis = [xmin, xmax, ymin, ymax, minZ, maxZ];
displayOption.caxis = [minZ, maxZ];



%% **************************************************************%

%-- sampling interval

dt =  0.35/signal.samplingFrequency*2;
Nt = round(signal.signalDuration/dt);
recordN = round(Nt/3)*2^(ref.maxRefine)-3;

%-- Source location

pSOR = xmin:dx_sensor:xmax;
%5IdS =  round(length(pSOR)/2); %test single emitter
pSORt = [pSOR(3),ymin];

%-- domain description

%load('rect_uni2.mat'); %uniform equilateral triangle
mesh.kABC = 3;                  % number of abosrbing edges
mesh.quadorder = 2;
expS   = strcat('y==',num2str(ymin));
expAB1 = strcat('x==', num2str(xmax));
expAB2 = strcat('y==',num2str(ymax));
expAB3 = strcat('x==',num2str(xmin));


%% Data Structure
%-- mesh & boundary process
%[mesh.node,mesh.elem] = squaremesh([xmin,xmax,ymin,ymax], signal.dx);
load('.\MeshData\rect_uni2.mat');
mesh.node = node;
mesh.elem = elem;

[mesh,dt,Nt,Ke,Kn] = meshTune(mesh,ref,Nt,dt,soundSpeedWater);

[mesh.elem2dof,edge,~] = dofP2(mesh.elem); % convert to P2 elem

% compute mesh parameters
mesh.nC = size(mesh.node,1);  
mesh.nE = size(mesh.elem,1); 
nEd = size(edge,1);
mesh.Ndof = mesh.nC + nEd;

% boundary data structure
mesh.solType = 'e';
bda = boundaryDataStructure(mesh,expS,expAB1,expAB2,expAB3);

%-- Stiffness matrix & mass matrix

disp(strcat('meshsize is',num2str(mesh.nC)));

% assemble stiffness matrices for global and boundary
[A,As,area] = stiffMatrixWaveP2(mesh,bda);

% assemble mass matrices for global and boundary
[M,Me,Ms] = massMatrixWaveP2(mesh,bda,area,Ke,Kn,dt);


% Initial zero conditions
u0 = zeros(mesh.Ndof,1);
u = zeros(mesh.Ndof,1);

% Initial absorbing boundary condition
phi = cell(mesh.kABC,1);
phid = phi;
phidd = phi;
for i = 1:mesh.kABC
    phi{i} = zeros(2*(2*bda.NS(i)-1),1);
    phid{i} = zeros(2*(2*bda.NS(i)-1),1);
    phidd{i} = zeros(2*(2*bda.NS(i)-1),1);
end

us = cell(mesh.kABC,1);
ut = cell(mesh.kABC,1);
Kt = cell(mesh.kABC,1);


%% Forward solver %

[lambda,weight] = quadpts1(4);
pxy1 = lambda(1,1)*mesh.node(bda.EdgeSOR(1:end-1),:) ...
     + lambda(1,2)*mesh.node(bda.EdgeSOR(2:end),:);
pxy2 = lambda(2,1)*mesh.node(bda.EdgeSOR(1:end-1),:) ...
     + lambda(2,2)*mesh.node(bda.EdgeSOR(2:end),:);
pxy3 = lambda(3,1)*mesh.node(bda.EdgeSOR(1:end-1),:) ...
     + lambda(3,2)*mesh.node(bda.EdgeSOR(2:end),:);
psi = zeros(3,3);
psi(:,1) = (2*lambda(:,1)-1).*lambda(:,1);
psi(:,2) = (2*lambda(:,2)-1).*lambda(:,2);
psi(:,3) = 4*lambda(:,1).*lambda(:,2);
bIn = zeros(length(bda.EdgeSOR)-1,3);
bdSOR2node = bda.SORid + mesh.nC;

for n=1:Nt  
    
    uIn1 = pde.Nsource(pxy1,pSORt,n*dt,signal);
    uIn2 = pde.Nsource(pxy2,pSORt,n*dt,signal);
    uIn3 = pde.Nsource(pxy3,pSORt,n*dt,signal);
          
  
    bIn(:,1) = (weight(1)*psi(1,1)*uIn1 + weight(2)*psi(2,1)*uIn2...
             + weight(3)*psi(3,1)*uIn3).*bda.elSOR;
    bIn(:,2) = (weight(1)*psi(1,2)*uIn1 + weight(2)*psi(2,2)*uIn2...
             + weight(3)*psi(3,2)*uIn3).*bda.elSOR;
    bIn(:,3) = (weight(1)*psi(1,3)*uIn1 + weight(2)*psi(2,3)*uIn2...
             + weight(3)*psi(3,3)*uIn3).*bda.elSOR;

    b = accumarray([bda.EdgeSOR(1:end-1);bda.EdgeSOR(2:end);bdSOR2node],bIn(:),[mesh.Ndof 1]); 
     
    % absorbing boundary condition
%     for i = 1:mesh.kABC
%         us{i} = [u0(bda.EdgeABC{i});u0(bda.ABCid{i}+ mesh.nC)];
%         ut{i} =[(u(bda.EdgeABC{i}(1))- u0(bda.EdgeABC{i}(1)))/dt,...
%                 (u(bda.EdgeABC{i}(bda.NS(i)))- u0(bda.EdgeABC{i}(bda.NS(i))))/dt];
%         Kt{i} = [sqrt(Kn(bda.EdgeABC{i}(1))), sqrt(Kn(bda.EdgeABC{i}(bda.NS(i))))];
%     end
% 
%     [phi,phid,phidd,Aux] =  WaveABCP2(us,ut,phi,phid,phidd,Ms,As,...
%                                       bda.NS,bda.bd2node,Kt,...
%                                       dt,mesh.Ndof,mesh.kABC);
     Aux = zeros(mesh.Ndof,1);
    % wave forward solver
    u1 = (M+0.5*Me)\((b-A*u) + M *(2*u-u0) + 0.5*Me*u0 + Me*Aux*dt);
    u0 = u;
    u = u1;
    
    if  displayF 
        ps = figure(1);
        showsolution(mesh.node,mesh.elem,u(1:mesh.nC));
        axis(displayOption.axis);
        colorbar
        caxis(displayOption.caxis);
        view(displayOption.viewAngle)
        pause(0.01)
        
    end
    
     if n==570*2
        ps = figure(1);
        showsolution(mesh.node,mesh.elem,u(1:mesh.nC));
        axis(displayOption.axis);
        colorbar
        caxis(displayOption.caxis);
        view(displayOption.viewAngle)
        pause(0.01)
        %saveas(ps,'p22_fem1.eps','epsc2');
      end
  
end