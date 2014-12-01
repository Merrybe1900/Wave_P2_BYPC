function signalFiltered = PerformLowpassFiltering(signal, samplingFrequency, lowpassFrequencyLimit, lowpassFilterType)
%----------------------------------------------------------
% Copyright 2011 SonoView Acoustic Sensing Technologies
%
% Author(s): Ivana Arsic, Ivana Jovanovic
% Description: This function performs lowpass filtering 
% such the all the frequencies higher than the lowpassFrequencyLimit
% are filtered out. This function uses built-in Matlab function filtfilt()
% that performs zero-phase filtering.  
%----------------------------------------------------------
% Filter parameters
frequencyAtPassBandBorder = lowpassFrequencyLimit;
frequencyAtStopBandBorder = 1.15 * lowpassFrequencyLimit;
attenuationInPassBand = 1;
attenuationInStopBand = 70;

% Desing lowpass filter.
h = fdesign.lowpass(frequencyAtPassBandBorder, frequencyAtStopBandBorder, attenuationInPassBand, attenuationInStopBand, samplingFrequency); 
switch lowpassFilterType 
    case {'butter','cheby1','cheby2','ellip'}
        d = design(h,lowpassFilterType);
        [b,a] = sos2tf(d.sosMatrix, d.ScaleValues);  
    case 'equiripple'
        d = design(h,lowpassFilterType);
        b = d.Numerator;
        a = 1;
end

% Apply filtering.
signalFiltered = filtfilt(b, a, signal); 