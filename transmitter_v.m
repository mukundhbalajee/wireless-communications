%% Setup
clear; clc; close all;

sync_len = 1; % microseconds
fs = 200; %MHz
upsampling_rate = 20;

% freq_sync = ones(1, fs*sync_len/upsampling_rate); %2 us of frame sync for the receiver to get PLL lock
% t_sync = floor(2.*rand(1, fs*sync_len/upsampling_rate)).* 2 -1; %Pseudo-random BPSK
% frame_sync = floor(2.*rand(1, fs*sync_len/upsampling_rate)).* 2 -1; %Pseudo-random BPSK
load freq_sync;
load t_sync.mat;
load frame_sync.mat;

preamble = horzcat(freq_sync, t_sync, frame_sync);
impulse = horzcat(zeros(1, 100), 1, zeros(1, 100));
sine = horzcat(zeros(1, 100), sin(2*pi*5e6.*linspace(0, 50e-6, 500)), zeros(1, 100));
square = horzcat(zeros(1, 100), ones(1, 100), zeros(1, 100));
test_sig = horzcat(preamble, impulse, zeros(1,40), sine, zeros(1, 40), square);
test_sig = upsample(test_sig, upsampling_rate);

data = imread('images/shannon1440.bmp');
dims = size(data);
bits = reshape(data, dims(1)*dims(2), 1);
bits = bits .* 2 -1; %BPSK
bits = horzcat(preamble, bits');
bits_up = upsample(bits, upsampling_rate);
pulse = rcosdesign(0.5, 7, upsampling_rate, 'sqrt');

transmitsignal = conv(test_sig, pulse);
save("transmitsignal.mat", "transmitsignal");
time = length(transmitsignal)/200; %In microseconds

%% Transmit

xt = transmitsignal;
len = length(xt);

figure(1)
LargeFigure(gcf, 0.15); % Make figure large
clf

% Plot time domain signal
subplot(2,1,1);
hold on;
plot([0:len-1]/200e6*1e6, real(xt))
plot([0:len-1]/200e6*1e6, imag(xt))
hold off;
xlabel('t in microseconds')
ylabel('x(t)')
legend("real", "imag");
title('xbase(t)')
axis tight

% Plot frequency domain signal
subplot(2,1,2);
plot([-len/2+1:len/2]*200/len, 20*log10(abs(fftshift(1/sqrt(len)*fft(xt)))))
xlabel('DTFT frequency f in MHz')
ylabel('abs(X(f)) in dB')
axis tight
xline(11, 'red');
xline(-11, 'red');
xline(15, 'black');
xline(-15, 'black');
yline(-20, 'red');
yline(-40, 'black');

%% Received
load receivedsignal.mat
xt = receivedsignal;
len = length(xt);

figure(2)
LargeFigure(gcf, 0.15); % Make figure large
clf

% Plot time domain signal
subplot(2,1,1);
hold on;
plot([0:len-1]/200e6*1e6, real(xt))
plot([0:len-1]/200e6*1e6, imag(xt))
hold off;
title('ybase(t)')
xlabel('t in microseconds')
ylabel('x(t)')
legend("real", "imag");
axis tight

% Plot frequency domain signal
subplot(2,1,2);
plot([-len/2+1:len/2]*200/len, 20*log10(abs(fftshift(1/sqrt(len)*fft(xt)))))
xlabel('DTFT frequency f in MHz')
ylabel('abs(X(f)) in dB')
xline(11, 'red');
xline(-11, 'red');
xline(15, 'black');
xline(-15, 'black');
yline(-20, 'red');
yline(-40, 'black');
axis tight
