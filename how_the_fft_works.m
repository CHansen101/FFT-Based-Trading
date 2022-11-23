%% How the FFT Works
clear; clc;
f = 60; % 60 Hz 
T = 1/f; % Period of our signal
Fs = 1000;% % 1k samples/second
Ts = 1/Fs; % Time difference between samples
t = 0:Ts:4*T; % time vector
x = cos(2*pi*f*t); % Sinusoid with frequency of 60 Hz
subplot(2,1,1); plot(t,x); % Plot a few periods of our signal 
title(sprintf('Signal with frequency of %i Hz',f));
xlabel('time (seconds)'); ylabel('Amplitude');
% As the time vector increases, the closer our fft becomes to the FT of our
% signal

N = 1024; % A higher value for N will produce a finer resolution 
% Note that the fft is most efficient when N is a power of 2.
X = abs(fft(x,N)); X = 2*X(1,1:N/2); % Compute the single sided fft
f_w = linspace(0,0.5,N/2)*Fs; % Generate the frequency axis 
subplot(2,1,2); plot(f_w,X) % Plot the fourier transform
title('FFT of Signal'); xlabel('Hz'); ylabel('|X(f)|');