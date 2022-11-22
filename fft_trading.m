% FFT Algorithm Based Trading
% author: Christopher Hansen
% date: 11/22/22
% This program is currently under development. The goal is to design a
% trading algorithm that will make trades based on the Fourier Transform of
% certain stock data. 
%% 
clear; clc;

% === Obtaining the historical data ===
% Technically, the provided data has a resolution of 1 day, however, I want
% to develop the algorithm for a finer resolution, so I am going to assume
% the resolution is 1 minute. Therefore, the "sampling frequency" is 60s or
% 1 minute. 
hist_data = readtable("qqq.csv"); % grab historical data
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
%% Developing Foundation for FFT computing on a large scale
N = 1024; % Number of points in the fft
Fs = 1/60; % sampling frequency in Hz
Ts = 1/Fs; % time difference between samples (seconds)
t = 0:Ts:dayToSec(5); % time vector for 5 days
f = 1/dayToSec(2); % Freqency of 1 period every 2 days (1/2*86400 seconds = Hz)
data_stream = cos(2*pi*f*t); % dummy data stream
H = fft(data_stream,N);
f_w = linspace(0,2*pi,N); % radians/second
%f_w = Fs*f_w/N; % convert to Hz

% Plot the dummy data
subplot(2,1,1); plot(secToDay(t),data_stream);
title('(Dummy) Stock Data Price - Known Period of 2 days');
xlabel('Days'); ylabel('Price');

% Plot the Fourier Transform of dummy data
subplot(2,1,2); plot(f_w,abs(H));
title('Fourier Transform of Stock Data Price');
xlabel('\omega (radians/second)'); ylabel('Amplitude');
%%
for overall_time = 1:length(hist_data(:,1)) % loop through the whole data set
    % Once an hour, compute the fft
    

end
