% FFT Algorithm Based Trading
% author: Christopher Hansen
% date: 11/22/22
% This program is currently under development. The goal is to design a
% trading algorithm that will make trades based on the Fourier Transform of
% certain stock data. 

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
%% Developing Foundation for FFT computing on a large scale
clear; clc;
N = 2^16; % Number of points in the fft
Fs = 1/60; % sampling frequency in Hz
Ts = 1/Fs; % time difference between samples (seconds)
t = 0:Ts:dayToSec(5); % time vector for 5 days
f = 1/dayToSec(2); % Freqency of 1 period every half day
data_stream = cos(2*pi*f*t); % dummy data stream
H = abs(fft(data_stream,N)); H = 2*H(1,1:N/2);
f_w = linspace(0,0.5,N/2)*Fs; % Hz (1/s)

% Plot the dummy data
subplot(2,1,1); plot(secToDay(t),data_stream);
title(sprintf('(Dummy) Stock Data Price - Known Period of %2.2f days',secToDay(1/f)));
xlabel('Days'); ylabel('Price');

% Plot the Fourier Transform of dummy data
subplot(2,1,2); plot(1./secToDay(1./f_w),H);
title('Fourier Transform of Stock Data Price');
xlabel('Periods per Days'); ylabel('Amplitude');

%%  === Obtaining the historical data ===
% Technically, the provided data has a resolution of 1 day, however, I want
% to develop the algorithm for a finer resolution, so I am going to assume
% the resolution is 1 minute. Therefore, the "sampling frequency" is 60s or
% 1 minute. 
clear; clc;
hist_data = readmatrix("qqq.csv"); % grab historical data
hist_data = flip(hist_data);
% below are the indexes of the columns in the hist_data variable
% associated with the specific data
date_ind = 1;   % index for the dates of prices in hist_data
close_ind = 2;  % index for the closing value of prices in hist_data
volume_ind = 3; % index for the volume of trades in hist_data
open_ind = 4;   % index for the opening value of prices in hist_data

% The algorithm will need to run in pseudo-real time, so we cannot take the
% fft over the whole data set. 
% 1) Let's assume we are going to start at the forefront of the datastream
% 2) We are not going to continuously compute the fft, as the fund will
% have a capital of less than $25,000, so we will be limited on how many
% trades we can make anyway. 
% 3) My first approach is to compute the fft every hour, or every (60) 1
% minute intervals. 
% 4) We will start with just getting data, appending it to our growing
% historical data set, and performing analysis on that data set once an
% hour. 

N = 2^16; % Number of points in the fft
Ts = dayToSec(1); % time difference between samples (seconds)
Fs = 1/Ts; % sampling frequency in Hz
t = 1:Ts:length(hist_data(:,1))*Ts;
subplot(2,1,1); plot(t,hist_data(:,close_ind));
title(sprintf('Stock Price of QQQ over Time (resolution of %2.2f day(s))',secToDay(Ts)));

