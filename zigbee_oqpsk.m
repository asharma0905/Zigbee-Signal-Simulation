% This function convert text into binary and return two variable
% One in vector -binV other in string binS
% Input: 
%    text - Class char/string e.g text = 'Hello World'
% Output:
%    binV - Binary vector of class double
%    binS - Binary Strin  of class char/string

function [binV, binS] = text2bin(text)
binS = dec2bin(text, 16); % dec2bin() returns binary string vectors.
binS = binS'; % Number of characters = number oof UTF-16 binary encodings.
binS = binS(:)'; % Forms the row from the binary string
binV = binS-48; % Converted to numbers from characters.
end


psdu = randi([0,1], 16, 1);
[binV, binS] = text2bin('Quick brown fox jumps over the lazy dog 1234567890 ~`!@#$%^&*()_-+={}[]|\;’:”<>,./?');
% The encoding used for Zigbee protocol depends on Application layer.
% UTF-16 used here.

% The row vector converted into column vector by taking transpose. 
binV = binV';

% The samples per chip parameter defined for defining the sampling rate to
% be used.
spc = 8;

% lrwpanOQPSKconfig is the inbuilt function of MATLAB to generate the
% configuration object for OQPSK modulation 802.15.4 transmitter PHY layer.
% It accepts the center frequency (Band), samplesPerChip, PSDU length and 
% the sampling rate of the output waveform as the parameters.

oqpsk_cfg = lrwpanOQPSKConfig(Band=2450, SamplesPerChip=spc);

% lrwpanWaveformGenerator is the inbuilt function of MATLAB which takes the
% data payload and the OQPSK config object as parameters to generate the 
% output modulated passband waveform.

waveOQPSK = lrwpanWaveformGenerator(psdu, oqpsk_cfg);

% Using DigitalUpConverter to upconvert the baseband signal inorder to 
% observe the frequency peaks at the center frequency. It accepts parameters
% center frequency, bandwidth, interpolation factor and sample rate amongst
% others. Sample Rate * interpolation factor should be equal to or greater 
% than Center frequency.

upConv = dsp.DigitalUpConverter("CenterFrequency", 2.45*1e9, ...
    "Bandwidth", 2*1e6, "InterpolationFactor", 2000, "SampleRate", 8*1e6 ...
    );

xUp = upConv(waveOQPSK);

% Plot of the imaginary and real part of the generated output waveform
% which shows the frequency shifted passband waveform in the time domain.
% Plotting the PSD vs frequency plot for the amplitude-time values after 
% upconverting the baseband signal generated by the inbuilt function in
% MATLAB. This served as a reference for verifying the plot generated by the
% simulated baseband signal values generated via the C++ code.

%plot(waveOQPSK)
%plot(real(waveOQPSK))
sz = size(xUp);
sz = sz(1,1);
fs = 16*1e6;
Y = fft((xUp));
%plot(20 * log10(abs(Y)));

% % hamming() function used to form the windows for the frequency values on 
% the x axis.
X = hamming(floor(length(xUp)/10));
% plot(fs/sz*(-sz/2:sz/2-1), 20 * log10(abs(fftshift(Y))));

% pwelch() used to plot the Power Spectrum Density vs the frequency values 
% on the x axis. It takes the Y values, X values as vectors, the sampling
% rate*interpolation factor for determining the range of x values. It also
% accepts a string as the type of the spectrum to be plotted.
pwelch(xUp, X, [], [], 1000*fs, 'centered')
% axis on

%xlabel('real(waveOQPSK)')
%ylabel('imag(waveOQPSK)')

% Eye diagram plot shows the amplitude time plot of the inphase and
% quadrature components of the output OQPSK passband modulated signal which 
% takes the parameter of the snapshot of the waveOQPSK to be plotted in the 
% eye diagram and number of samples per chip.

%eyediagram(waveOQPSK(1+spc:end-spc), 2*spc);