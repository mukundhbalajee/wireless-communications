%% Setup
clear; clc; close all;
load receivedsignal.mat
data = imread('images/shannon1440.bmp');
dims = size(data);
bits = reshape(data, dims(1)*dims(2), 1);
bits = bits .* 2 -1; %BPSK
bits_up = upsample(bits, 20);
pulse = rcosdesign(0.5, 7, 20, 'sqrt');
transmitsignal = conv(bits_up, pulse);
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
