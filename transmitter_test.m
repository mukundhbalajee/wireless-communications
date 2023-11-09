% Basic BPSK Transmitter Design

% Load the image and convert it to a binary sequence
img = imread('images/shannon1440.bmp');
message = reshape(img, [], 1);

% Values are either '0' or '1' -> multiple values by 2 
% and subtract 1 to make them '-1' and '1', to create message vector
xk = (message * 2) - 1;

% Parameters
F_samp = 200e6; % Sampling rate, 200 MHz
T = 1/(22e6); % Sampling period
%floor to avoid going over the limits
fs = floor(F_samp*T);% Over-sampling factor (Sampling frequency/symbol rate).
N = 51; % Length of filter in symbol periods. Default is 51
f_c = 11e6; % Carrier frequency, 11 MHz based on mask

Ns = floor(N*fs); % Number of filter samples

% Choose sinc pulse
% '1/fs' simply serves as 'delta' to approximate integral as sum
pt = sinc([-floor(Ns/2):Ns-floor(Ns/2)-1]/fs); 
pt=transpose(pt)/norm(pt)/sqrt(1/(fs)); 

% Create baseband signal
xk_up = upsample(xk,fs);

% Create the baseband signal (complex-valued)
x_I = real(xk_up); % In-phase component
x_Q = imag(xk_up); % Quadrature component

x_I = conv(x_I,pt);
x_Q = conv(x_Q,pt);

x_base = x_I + 1i*x_Q;

% Ensure signal magnitude is within [-1, 1]
if max(abs(x_base)) > 1
    x_base = x_base ./ max(abs(x_base));
end

% Save the signal to a .mat file
transmitsignal = x_base;
save('transmitsignal.mat', 'transmitsignal');

transmitted_signal = [1:length(transmitsignal)]/F_samp*1e6;
len = length(transmitted_signal);

% Plotting
figure;

% Time Domain
subplot(2,1,1);
plot(transmitted_signal, real(transmitsignal),'b')
hold on
plot(transmitted_signal, imag(transmitsignal),'r')
legend('real','imag')
ylabel('x(t)')
xlabel('Time in microseconds')

% Frequency Domain
subplot(2,1,2);
f_mhz = (-len/2+1:len/2) / len * F_samp / 1e6;
magnitude = 20*log10(abs(fftshift(1/sqrt(len)*fft(transmitsignal))));
plot(f_mhz, magnitude);
xlabel('Frequency (MHz)');
ylabel('|X(f)| (dB)');
title('Baseband Signal in Frequency Domain');
