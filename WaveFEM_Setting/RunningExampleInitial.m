function [signal,ref,displayOption,mesh,source,Nt,dt,maxCount,expS,...
                         expAB1,expAB2,expAB3,sf_Rec,sf_Err,sf_Mesure2,...
                         ini_file, second_ini_file] = RunningExampleInitial
%************************ User setting function **************************%
% this function defines all the user setting for the whole reconstruction
% procedure.
%*************************************************************************%
% Author: Yun Bai
% Date: 22.08.2014
%*************************************************************************%
soundSpeedWater = 1500; % [m/s]
pde = SignalSetup;


%**************************** Test setting ******************************%
%--- Computing problem ----%
isReflectingTest = false;   %in this case expAB2 is reflecting boundary. 
source = @pde.DNsource;    % the type of source signal is used
signal.transform = false;  %if Hilbert transform is applied
signal.p = 0.6;            % specify the parameter p in Hilbert transform
maxCount = [250,150];      % if transform is applied indicate [maxiter1,maxiter2]
                           % otherwise maxiter.

ref.ishomo = 0;            % if the reference phantom is a test phantom
                           % if ref.ishomo = 0, indicate 
                           % ref.testType = 6 (big phantom) or 7 (small phantom)
ref.testType = 6;          % if the ref.ishomo = 1 then need to indicate
                           % ref.testType = 1,2,3,4,5
                           
ref.isReceivedSignal = 1;          % if the forward problem needs to compute                            

%-- Mesh setting ---%
ref.maxRefine = 1;                 % refinement times
ref.isUniform = 1;                 % if it is uniform refinement
% ref.Nref = 2;                    % if refinement is nonuniform, then need to indicate the times
                                   % of nonuniform refinement
dx = 2.0000e-03;         % initial coarse meshsize [m]

%-- Sensor spacial position ----%
signal.dx = 2.0000e-03;  % spacial size of the signal  [m]
dx_sensor = 3*dx;        % the distance between two emitters
                         % (we assume all the nodes on the emitting BC
                         %  receive signals but a few of them emit)
        
%--- save file Name ---%
sf_Rec = '.\WaveRecData\reconstruction_P2_PE_';
sf_Err = '.\WaveRecData\reconstruction_P2_PE_CostFunction.mat';

%*************************************************************************%

%--- Domain description ---%
if ref.testType == 7
    xmin = -0.05;
    xmax = 0.05;
    ymin = -0.005;  % emitting b.c.
    ymax = 0.035;
elseif ref.testType == 6
    xmin = -0.08;
    xmax = 0.08;
    ymin = -0.022;  % emitting b.c.
    if ~isReflectingTest 
        ymax = 0.042;
    else % reflection test is only on bigger phantom
        ymax = 0.036 + 0.006*rand(1); 
    end
end


%--- Main Signal description ---%
signal.signalDuration = 180e-6; % [s]
signal.samplingFrequency = 5e6; % [Hz]
signal.sigtype = 'exp';  % signal expression 'exp' or 'sinc'

signal.t0 = 2e-6;        % initial delay [s]
signal.am = 100;          % signal amplitude 
signal.fc = 150e3;       % central frequency [Hz]
signal.bw = 150e3;       % bandwidth [Hz]

signal.sensorNum = length(signal.fc);  % number of sensor types
signal.pSOR = (xmin+ dx:dx_sensor:xmax-9*dx)'; % emitters' locations

signal.numS = length(signal.pSOR);              % number of emitters


% if transform is used, the transformed received signal is used 
% received signal file
signal.sf_Mesure = '.\WaveMesureData\receivedSignal.mat';

dt =  0.3/signal.samplingFrequency*2;           % time step size
Nt = round(signal.signalDuration/dt);           % number of total time steps
disp(strcat('total number of sensor is ', num2str(signal.numS)));


%-- Reference forward setting (or for solver test) ---%
load('.\WaveMeshData\rect_uni2.mat');  % initial mesh for reference forward 
% for the forward problem, the computing domain is fixed to either smaller
% phantom or bigger one
%[node,elem] = squaremesh([xmin,xmax,ymin,ymax], dx);
ref.node = node;                   % reference mesh node cordinates
ref.elem = elem;                   % reference element node indices

ref.quadorder = 3.5;               % 3.5: masslumpting p2-FEM
                                   % 2: classcial p2-FEM
ref.sf_Mesh = '.\WaveMeshData\';
if isReflectingTest 
   ref.rib_opt = 0;               % then the domain should not contain ribs
   ref.kABC = 2;                      % number of absorbing edges
else
   ref.rib_opt = 1;               % then the domain should not contain ribs
   ref.kABC = 3;                      % number of absorbing edges
end

%--- Display option ---%
displayOption.displayF = 0;         % Forward solution display
displayOption.displayB = 0;         % Backward solution display
displayOption.viewAngle = [0,-90];
maxZ = 0.08;  minZ = -0.03;         % Display range
displayOption.axis = [xmin, xmax, ymin, ymax, minZ, maxZ];
displayOption.caxis = [minZ, maxZ];


%--- reconstruction domain description ----%
expS   = strcat('y==',num2str(ymin));
expAB1 = strcat('x==', num2str(xmax));
expAB2 = strcat('y==',num2str(ymax));
expAB3 = strcat('x==',num2str(xmin));


%--- Parameters for reconstruction ---%
mesh.kABC = ref.kABC;                      % number of absorbing edges
mesh.maxIterL = 3;       % max iterations for line search
mesh.quadorder = ref.quadorder;    % 3.5: masslumpting p2-FEM
                         % 2: classical p2-FEM
mesh.K0 = 1/soundSpeedWater;
mesh.ymin = ymin;
mesh.rangeY = ymax-ymin;

% initial mesh for reconstrution, normally different from reference mesh 
[mesh.node,mesh.elem] = squaremesh([xmin,xmax,ymin,ymax], dx); % right triangle uniform mesh

% file names
sf_Mesure2 = '.\WaveMesureData\reflectionSignal.mat';
if signal.transform
   ini_file = [];
   second_ini_file = strcat(sf_Rec,num2str(maxCount(1)),'.mat');
else
   ini_file = [];
   second_ini_file = [];
end
end