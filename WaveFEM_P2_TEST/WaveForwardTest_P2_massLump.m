function WaveForwardTest_P2_massLump
soundSpeedWater = 1500; % [m/s]
K0 = 1/soundSpeedWater; % constant 1/c0
pde = SignalSetup;

%--- signal description ---%
signal.signalDuration = 120e-6; % [s]
signal.samplingFrequency = 5e6; % [Hz]
signal.fc = 150e3; % [Hz]
signal.bw = 150e3;       % [Hz]
signal.sigtype = 'exp';      % 'exp' or 'sinc'
signal.t0 = 2e-6;        % [s]
signal.am = 120;
signal.dx = 2.0000e-03; 
%pde = SignalSetup;

%-- Computational domain description ---%
xmin = -0.08;
xmax = 0.08;
ymin = -0.022;
ymax = 0.042;
        % initial coarse meshsize
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
pSORt = [pSOR(3),ymin]

%-- domain description
%load('rect_uni2.mat'); %uniform equilateral triangle
mesh.kABC = 3;                  % number of abosrbing edges
mesh.quadorder = 3.5;
expS   = strcat('y==',num2str(ymin));
expAB1 = strcat('x==', num2str(xmax));
expAB2 = strcat('y==',num2str(ymax));
expAB3 = strcat('x==',num2str(xmin));


%% data structure

%-- mesh & boundary process
%[mesh.node,mesh.elem] = squaremesh([xmin,xmax,ymin,ymax], dx);
load('.\MeshData\rect_uni2.mat');
mesh.node = node;
mesh.elem = elem;
mesh.solType = 'e';
[mesh,dt,Nt,Ke,Kn] = meshTune(mesh,ref,Nt,dt,soundSpeedWater);

[mesh.elem2dof,edge,~] = dofP2(mesh.elem); % convert to P2

% compute mesh parameters
mesh.nC = size(mesh.node,1);  
mesh.nE = size(mesh.elem,1); 
nEd = size(edge,1);
mesh.elem2dof = [mesh.elem2dof,(1+nEd+mesh.nC:mesh.nC + nEd + mesh.nE)'];
mesh.Ndof = mesh.nC + nEd + mesh.nE;%-- Boundary process

% boundar data structure
bda = boundaryDataStructure(mesh,expS,expAB1,expAB2,expAB3);



%-- Stiffness matrices & mass matrices
disp(strcat('meshsize is',num2str(mesh.nC)));

[A,As,area] = stiffMatrixWaveP2(mesh,bda,Kn);
[M,Me,ML] = massMatrixWaveP2(mesh,bda,area,Ke,Kn,dt);

% Initial zeros conditions
u0 = zeros(mesh.Ndof,1);
u1 = zeros(mesh.Ndof,1);


%%  Forward solver %

lambda = [1,0,0.5;0,1,0.5]';
weight = [1/6,1/6,2/3];
pxy1 = mesh.node(bda.EdgeSOR(1:end-1),:);
pxy2 = mesh.node(bda.EdgeSOR(2:end),:);
pxy3 = 0.5*(pxy1 + pxy2);   

psi = zeros(3,3);
psi(:,1) = (2*lambda(:,1)-1).*lambda(:,1);
psi(:,2) = (2*lambda(:,2)-1).*lambda(:,2);
psi(:,3) = 4*lambda(:,1).*lambda(:,2);
bIn = zeros(length(bda.EdgeSOR)-1,3);
bdSOR2node = bda.SORid + mesh.nC;

for n = 2:3 % first order ABC to initialize
    uIn1 = pde.Nsource(pxy1,pSORt,n*dt,signal);
    uIn2 = pde.Nsource(pxy2,pSORt,n*dt,signal);
    uIn3 = pde.Nsource(pxy3,pSORt,n*dt,signal);
          
  
    bIn(:,1) = (weight(1)*psi(1,1)*uIn1 + weight(2)*psi(2,1)*uIn2...
             + weight(3)*psi(3,1)*uIn3).*bda.elSOR;
    bIn(:,2) = (weight(1)*psi(1,2)*uIn1 + weight(2)*psi(2,2)*uIn2...
             + weight(3)*psi(3,2)*uIn3).*bda.elSOR;
    bIn(:,3) = (weight(1)*psi(1,3)*uIn1 + weight(2)*psi(2,3)*uIn2...
             + weight(3)*psi(3,3)*uIn3).*bda.elSOR;
         
    b = accumarray([bda.EdgeSOR(1:end-1);bda.EdgeSOR(2:end);bdSOR2node],...
                    bIn(:),[mesh.Ndof 1]); 
     
    if n ==2
       u2 =  (b  + (2*M-A)*u1 + (0.5*Me - M)*u0).*ML;
    elseif n == 3
       u3 =  (b  + (2*M-A)*u2 + (0.5*Me - M)*u1).*ML;
    end
     
end

for n=4:Nt-1  % secdon order mass lumption 

    uIn1 = (pde.Nsource(pxy1,pSORt,(n+1)*dt,signal)...
         -  pde.Nsource(pxy1,pSORt,(n-1)*dt,signal));
    uIn2 = (pde.Nsource(pxy2,pSORt,n*dt,signal)...
         - pde.Nsource(pxy2,pSORt,(n-1)*dt,signal));
    uIn3 = (pde.Nsource(pxy3,pSORt,n*dt,signal)...
         -  pde.Nsource(pxy3,pSORt,(n-1)*dt,signal));

    bIn(:,1) = (weight(1)*psi(1,1)*uIn1 + weight(2)*psi(2,1)*uIn2...
             + weight(3)*psi(3,1)*uIn3).*bda.elSOR;
    bIn(:,2) = (weight(1)*psi(1,2)*uIn1 + weight(2)*psi(2,2)*uIn2...
             + weight(3)*psi(3,2)*uIn3).*bda.elSOR;
    bIn(:,3) = (weight(1)*psi(1,3)*uIn1 + weight(2)*psi(2,3)*uIn2...
             + weight(3)*psi(3,3)*uIn3).*bda.elSOR;


     b = accumarray([bda.EdgeSOR(1:end-1);bda.EdgeSOR(2:end);bdSOR2node],...
                    bIn(:),[mesh.Ndof 1]); 
    % wave forward solver
    u4 =  (b  + (2*M-A)*u3 + (Me - 2*dt*As)*u2 + (A-2*M)*u1 - (0.5*Me-M)*u0).*ML;
    u0 = u1;
    u1 = u2;
    u2 = u3;
    u3 = u4; 

    %     
    if  displayF 
        ps = figure(1);
        showsolution(mesh.node,mesh.elem,u4(1:mesh.nC))
        axis(displayOption.axis);
        colorbar
        caxis(displayOption.caxis);
             view(displayOption.viewAngle)
        pause(0.01)
        
        if n==570*2
           saveas(ps,'p2_fem1.eps','epsc2');
        end
    end
end