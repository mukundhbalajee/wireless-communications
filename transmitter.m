%% Setup
clear; clc; close all;

sync_len = 24; % microseconds
fs = 200; %MHz
upsampling_rate = 12;

new_preamble = 0;

if(new_preamble)
    disp("New Preamble")

    freq_sync = ones(1, fs*sync_len/upsampling_rate); %2 us of frame sync for the receiver to get PLL lock
    t_sync = floor(2.*rand(1, fs*sync_len/upsampling_rate)); %Pseudo-random
    frame_sync = floor(2.*rand(1, fs*sync_len/upsampling_rate)); %Pseudo-random
    save("frame_sync", "frame_sync")
    save("t_sync", "t_sync")
    save("freq_sync", "freq_sync")
else
    load freq_sync;
    load t_sync.mat;
    load frame_sync.mat;
end

preamble = horzcat(freq_sync, t_sync, frame_sync);

if(mod(length(frame_sync), 4)) error("NOT MOD 4"); end


%% Image to QAM

image = imread('images/shannon20520.bmp');
dims = size(image);
im = double(reshape(image, 2052, dims(1)*dims(2)/2052));

bits = horzcat(preamble, im(:,1)', frame_sync, im(:,2)', frame_sync, im(:,3)', ...
    frame_sync, im(:,4)', frame_sync, im(:,5)', frame_sync, im(:,6)', frame_sync, ...
    im(:,7)', frame_sync, im(:,8)', frame_sync, im(:,9)', frame_sync, im(:,10)');

bits = string(reshape(bits, 4, length(bits)/4));

qam_bits = strings(1, length(bits));
for i = 1:length(bits)
    qam_bits(i) = strjoin(bits(:,i),'');
end

qam_bits = bin2dec((qam_bits))+1;

d = sqrt(2)/3;
options = [1.5+1.5j, 0.5+1.5j, -1.5+1.5j, -0.5+1.5j
           1.5+0.5j, 0.5+0.5j, -1.5+0.5j, -0.5+0.5j
           1.5-1.5j, 0.5-1.5j, -1.5-1.5j, -0.5-1.5j
           1.5-0.5j, 0.5-0.5j, -1.5-0.5j, -0.5-0.5j].*d;

qam_points = options(qam_bits);

%% Transmission

bits_up = upsample(qam_points, upsampling_rate);
pulse = rcosdesign(0.3, 30, upsampling_rate, 'sqrt');
save('bits', 'bits');

transmitsignal = conv(bits_up, pulse);
save("transmitsignal.mat", "transmitsignal");
time = length(transmitsignal)/200; %In microseconds

%% Transmitted signal - x_base(t)

x_base_t = transmitsignal;
len = length(x_base_t);
t_microseconds=[0:len-1]/200e6*1e6;

figure(1)
LargeFigure(gcf, 0.15); % Make figure large
clf

% Plot time domain signal
subplot(2,1,1);
hold on;
plot(t_microseconds, real(x_base_t))
plot(t_microseconds, imag(x_base_t))
hold off;
xlabel('t in microseconds')
ylabel('x^{base}(t)')
legend("real", "imag");
title('Time Domain Plot for transmitted signal, x_{base}(t)')
axis tight

% Plot frequency domain signal
subplot(2,1,2);
plot([-len/2+1:len/2]*200/len, 20*log10(abs(fftshift(1/sqrt(len)*fft(x_base_t)))))
xlabel('DTFT frequency f in MHz')
ylabel('|X^{base}(f)| in dB')
title('Frequency Domain Plot for transmitted signal, x_{base}(t)')
xline(11, 'red');
xline(-11, 'red');
xline(15, 'black');
xline(-15, 'black');
yline(-20, 'red');
yline(-40, 'black');
axis tight

