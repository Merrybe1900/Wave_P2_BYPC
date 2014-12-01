% function signalTest
t0 = 2e-6;        % initial delay [s]
signal.am = 100;          % signal amplitude 
fc = 150e3;       % central frequency [Hz]
bw = 150e3;       % bandwidth [Hz]

signal.signalDuration = 180e-6; % [s]
signal.samplingFrequency = 5e6; % [Hz]
sigtype = 'exp';  % signal expression 'exp' or 'sinc'
dt =  0.3/signal.samplingFrequency*2;           % time step size
Nt = round(signal.signalDuration/dt);           % number of total time steps
dtS =(dt:dt:Nt*dt)';
s = generatesource(fc, bw, dtS, sigtype, t0);
sH = abs(hilbert(s));

Hil = 1./(dtS) *(1/pi);
Hilc = conv(Hil,s,'same');
figure;

subplot(3,1,1);
plot((dt:dt:Nt*dt)',s);
subplot(3,1,2);
plot((dt:dt:Nt*dt)',sH);
subplot(3,1,3);
plot((dt:dt:Nt*dt)',Hilc);
% end