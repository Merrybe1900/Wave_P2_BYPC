function pde = SignalSetup
% Input parameters
%   p:       source boundary coordinates (the one coordinate which is
%            changing -- e.g. x, if y is constant) [m]
%   t:       current time instant [s]
%   fc:      center frequency [Hz]
%   bw:      bandwidth [Hz]
%   sigtype: signal type ('sinc' or 'exp')
%   t0:      initial delay [s]
%   am:      signal amplitude

pde = struct('Nsource', @Neumannsource,...
             'ENsource',@ExtendNsource,...
             'DNsource',@DelayNsource);

%----- subfunctions -----------%
% point source
    function z = Neumannsource(p,pSOR,t,signal) % Point source

           spaceF = exp(-(p(:,1)-pSOR(1)).^2/(1.4*(signal.dx^2)));
           timeF = generatesource(signal.fc, signal.bw, t)*signal.am;
           z =  spaceF.*timeF;

           
    end

% extended source
    function z = ExtendNsource(p,pSOR,t,signal)  % space extended source
            timeF = generatesource(signal.fc, signal.bw, t)*signal.am;
            indicator = zeros(length(p),1);
            indicator(abs(p(:,1)-pSOR(:,1))<= 0.005) = 1;
            z =  indicator.*timeF;
    end

% delayed source
    function z = DelayNsource(p,pSOR,t,signal) % series delayed source
           z = zeros(length(p),1);
           for i = 1:8
               spaceF = exp(-(p(:,1)-pSOR(1)-(i-1)*signal.dx).^2/(1.4*(signal.dx^2)));
               % delayed angle 30 degress
               timeF = generatesource(signal.fc, signal.bw, t,'exp', ...
                    1.5/signal.bw + (i-1)*signal.dx/1500*2/sqrt(3))*signal.am;
               z =  z + spaceF.*timeF;
           end
     end
end