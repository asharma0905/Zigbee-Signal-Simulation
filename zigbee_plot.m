% readMatrix() used to read the data of csv into a matrix.
A = readmatrix('/Users/amartyasharma/Downloads/Zigbee_C++/test_spc.csv');
a = A(:,2);
%[minA, maxA] = bounds(A(:,2));
% a = a/maxA;
b = A(:,3);
%[minb, maxb] = bounds(A(:,3));
% b = b/maxb;
%Y = [a, b];
%Y = repmat(Y, [100,1]);
compl = complex(a,b);

% plot(A(:,1), A(:,2))

% size() used to get the length of vector A into a matrix [len(A), 1] 
sz = size(A);

% Accessing the elements of A.
sz = sz(1,1); 
disp(sz);

% Using DigitalUpConverter to upconvert the baseband signal inorder to 
% observe the frequency peaks at the center frequency. It accepts parameters
% center frequency, bandwidth, interpolation factor and sample rate amongst
% others. Sample Rate * interpolation factor should be equal to or greater 
% than Center frequency.

upConv = dsp.DigitalUpConverter("CenterFrequency", 2.45*1e9, ...
"Bandwidth", 2*1e6, "InterpolationFactor", 1000, "SampleRate", 16*1e6 ...
);
Y = upConv(compl);

% bounds() function used to return min and max values of Y as a vector.
%[minY,maxY] = bounds(Y);
% Plotting the columns of A as x and y.
% plot(A(:,1), A(:,2));

% Plotting the DFT of the complex amplitude values to observe upconverted
% signal in frequency domain.

% fft(X) where X is a column vector of values formed by the amplitude time
% plot.
L = sz;
% Y = fft(compl);
fs = 16*1e6;
% X = fs/L*(0:L-1);

% Y = downsample(Y, 1000);

% hamming() function used to form the windows for the frequency values on 
% the x axis. 
X = hamming(floor(length(compl)/10));

% pwelch() used to plot the Power Spectrum Density vs the frequency values 
% on the x axis. It takes the Y values, X values as vectors, the sampling
% rate*interpolation factor for determining the range of x values. It also
% accepts a string as the type of the spectrum to be plotted.
pwelch(Y,X,[],[],1000*fs, "centered");

% Code to generate the .wav file by writing the generated Y values after
% upconversion to the carrier frequency. Scaling performed to limit the
% values between -1 to 1 for double datatype values in the vector Y. 
filename = 'zigbee_bb.wav';
%scaled_Y = Y/maxY;
%audiowrite(filename,Y,fs)

% Plotting the real part of FFT.
% plot(fs/L*(-L/2:L/2-1), 20 * log10(abs(fftshift(Y)))); 

% Plotting the imaginary part of FFT.
% plot(X,abs(Y)); 