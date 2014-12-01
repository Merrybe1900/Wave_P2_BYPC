function showSignal

numSor = 1;
numRec = 8;


signal.signalDuration = 180e-6; % [s]
signal.samplingFrequency = 5e6; % [Hz]
dt =  0.3/signal.samplingFrequency*2;          % time step size
Nt = round(signal.signalDuration/dt)*2;           % number of total time steps

sf_Mesure = '.\WaveMesureData\receivedSignal.mat';
sf_Mesure2 = '.\WaveMesureData\reflectionSignal.mat';

load(sf_Mesure,'uSOR');

figure(1);
plot(dt:dt:Nt*dt,uSOR{numSor}(numRec,:),'k');
ylim([-0.01,0.01]);
load(sf_Mesure2,'uhS');

figure(2);
plot(dt:dt:Nt*dt,uhS{numSor}(numRec,:),'r');
ylim([-0.01,0.01]);


end