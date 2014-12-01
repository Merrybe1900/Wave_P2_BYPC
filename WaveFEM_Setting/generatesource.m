function s = generatesource(fc, bw, t, sigtype, t0)
% generatesource - generate a signal with a given center frequency,
% bandwidth, sampling frequency and shape
%
% Input parameters:
%   fc:      center frequency of the bandpass pulse [Hz]
%   bw:      20-dB bandwidth of the bandpass pulse [Hz]
%   t:       vector of time samples (sampling time instants) [s]
%   sigtype: string which determines the pulse shape
%            'exp' (default): exponential function modulated by a cosine
%            'sinc':          difference between two sincs (not
%                             recommended due to slow decays)
%   t0:      pulse delay [s] (default = 1.5 / bw)
%
% Outputs:
%   s:       generated signal
%
% Example:
%   s = generatesource(100000, 50000, 0:1e-6:800e-6, 'exp', 80e-6)
%   s = generatesource(100000, 50000, 0:1e-6:800e-6, 'sinc')
%   s = generatesource(100000, 50000, 0:1e-6:800e-6)


% Copyright 2013-2014 SonoView Acoustic Sensing Technologies.
%
% Author:  Mihailo Kolundzija
% Created: 26-Feb-2014

if nargin < 5
  t0 = 1.5 / bw; % pulse delay
end

if nargin < 4
  sigtype = 'exp';
end

if strcmp(sigtype, 'sinc') == 1
  f0 = fc-bw/2;
  f1 = fc+bw/2;
  s = 2*pi*sinc(2*f1*(t-t0)) - 2*pi*f0/f1*sinc(2*f0*(t-t0));
elseif strcmp(sigtype, 'exp') == 1
  tc = 1 / (2.22*bw);
  t0 = 3.6 * tc;
  tt = t - t0;
  s = exp(-(tt .* tt) / (2*tc.*tc)) .* cos(2*pi*fc.*tt);
else
  error(['Invalid signal type: ' sigtype]);
end
