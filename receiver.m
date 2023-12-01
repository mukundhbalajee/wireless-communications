clear; close all; clc;
load receivedsignal.mat;
load transmitsignal.mat;
load qam_pre.mat
load qam_fsync.mat

preamble = qam_preamble;
frame_sync = qam_frame_sync;

sync_len = 24; % microseconds
fs = 200; %MHz
upsampling_rate = 12;
pulse = rcosdesign(0.3, 30, upsampling_rate, 'sqrt');

%% Demodulate
wt = fliplr(pulse); %probably unnecessary
zt = conv(wt, receivedsignal);

%% Timing Synchronization

sig_to_find = preamble;  %PUT SIGNAL TO FIND HERE

xt_conj = conj(sig_to_find); 

tau = length(zt);
sums = zeros(1, tau);
p = length(xt_conj);

for i = 1:tau-p*upsampling_rate
    sums(i) = 0;
    for k = 0:p-1
        sums(i) = sums(i) + xt_conj(k+1)*zt(k*upsampling_rate+i);
    end
end

timing_offset = find(abs(sums) == max(abs(sums)));
zt_timing = zt(timing_offset+upsampling_rate:end);
phase = angle(max(sums));
zt_inphase_timing = zt_timing .* exp(-j*phase); % Remove phase offset
zt_inphase_timing_plot = zt_inphase_timing;

%% Frame Synchronization

sig_to_find = frame_sync;  %PUT SIGNAL TO FIND HERE
num_frames = 10;
xt_conj = conj(sig_to_find);

tau = length(zt_inphase_timing);
sums = zeros(1, tau);
p = length(xt_conj);

for i = 1:tau-p*upsampling_rate-1
    sums(i) = 0;
    for k = 1:p
        sums(i) = sums(i) + xt_conj(k)*zt_inphase_timing(k*upsampling_rate+i);
    end
end

[pks,inds] = findpeaks(abs(sums), 'MinPeakProminence', 5);
maxes = maxk(pks, num_frames);

locs = zeros(1, num_frames);

for i = 1:num_frames
    locs(i) = find(abs(sums) == maxes(i));
end

locs = sort(locs);

phases = zeros(1, num_frames);
starts = zeros(1, num_frames);
for i = 1:length(inds)
   phases(i) = angle(sums(locs(i)));
   starts(i) = locs(i) + upsampling_rate*(length(frame_sync)+1);
end

%% Equalization
received_preamble = zt_inphase_timing(1:upsampling_rate:length(preamble)*upsampling_rate);
channel_effect = abs(received_preamble)'./abs(preamble);
h0=mean(abs(channel_effect));
zt_inphase_timing = zt_inphase_timing./h0;

%% Sampling
% zk = zt_inphase_frame(1:upsampling_rate:length(zt_inphase_frame));
frame1 = zt_inphase_timing(starts(1): starts(2)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(1));
frame2 = zt_inphase_timing(starts(2): starts(3)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(2));
frame3 = zt_inphase_timing(starts(3): starts(4)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(3));
frame4 = zt_inphase_timing(starts(4): starts(5)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(4));
frame5 = zt_inphase_timing(starts(5): starts(6)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(5));
frame6 = zt_inphase_timing(starts(6): starts(7)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(6));
frame7 = zt_inphase_timing(starts(7): starts(8)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(7));
frame8 = zt_inphase_timing(starts(8): starts(9)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(8));
frame9 = zt_inphase_timing(starts(9): starts(10)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(9));
frame10 = zt_inphase_timing(starts(10): end-1-upsampling_rate*(length(frame_sync))).* exp(-j*phases(10));

zk1 = frame1(1:upsampling_rate:length(frame1));
zk2 = frame2(1:upsampling_rate:length(frame2));
zk3 = frame3(1:upsampling_rate:length(frame3));
zk4 = frame4(1:upsampling_rate:length(frame4));
zk5 = frame5(1:upsampling_rate:length(frame5));
zk6 = frame6(1:upsampling_rate:length(frame6));
zk7 = frame7(1:upsampling_rate:length(frame7));
zk8 = frame8(1:upsampling_rate:length(frame8));
zk9 = frame9(1:upsampling_rate:length(frame9));
zk10 = frame10(1:upsampling_rate:length(frame10));

zk = [zk1; zk2; zk3; zk4; zk5; zk6; zk7; zk8; zk9; zk10];
zk = transpose(zk);


% Guessing
d = sqrt(2)/3;
options = [1.5+1.5j, 0.5+1.5j, -1.5+1.5j, -0.5+1.5j
           1.5+0.5j, 0.5+0.5j, -1.5+0.5j, -0.5+0.5j
           1.5-1.5j, 0.5-1.5j, -1.5-1.5j, -0.5-1.5j
           1.5-0.5j, 0.5-0.5j, -1.5-0.5j, -0.5-0.5j].*d;

received_pts = zeros(1, length(zk)); %QAM

min_pt = -1;
min_distance = 10000;
for i = 1:length(zk)
    min_pt = -1;
    min_distance = 10000;
    for j = 1:16
        refpoint = [real(options(j)), imag(options(j))];
        point = [real(zk(i)), imag(zk(i))];
        distance = sqrt((refpoint(1)-point(1))^2 + (refpoint(2)-point(2))^2);
        if(distance < min_distance)
            min_distance = distance;
            min_pt = j;
        end
    end
    received_pts(i) = min_pt-1;
end

% %QAM to Bits
% rx_bits = dec2bin(received_pts);
% rx_bits = string(rx_bits);
% rx_bits = strjoin(rx_bits,'');
% rx_bits = strsplit(rx_bits);


%Error checking and BER calculation
image = imread('images/shannon20520.bmp');
dims = size(image);

bits = string(double(reshape(image, 4, dims(1)*dims(2)/4)));

qam_bits = strings(1, length(bits));
for i = 1:length(bits)
    qam_bits(i) = strjoin(bits(:,i),'');
end

qam_bits = bin2dec((qam_bits));
N = length(qam_bits); %length(qam_bits);
BER= length(find(received_pts(1:N) ~= qam_bits(1:N)))/N*100

%% Plotting

% figure; imshow(data); title("original")
% 
% received_image = [rx_ipts1(1:400); rx_ipts2(1:400); rx_ipts3(1:400); rx_ipts4(1:240)];
% received_image = reshape(received_image, 45, 32);
% figure; imshow(received_image); 


% %% Transmitted signal - x_base(t)
% 
% x_base_t = transmitsignal;
% len = length(x_base_t);
% t_microseconds=[0:len-1]/200e6*1e6;
% 
% figure(1)
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% 
% % Plot time domain signal
% subplot(2,1,1);
% hold on;
% plot(t_microseconds, real(x_base_t))
% plot(t_microseconds, imag(x_base_t))
% hold off;
% xlabel('t in microseconds')
% ylabel('x^{base}(t)')
% legend("real", "imag");
% title('Time Domain Plot for transmitted signal, x_{base}(t)')
% axis tight
% 
% % Plot frequency domain signal
% subplot(2,1,2);
% plot([-len/2+1:len/2]*200/len, 20*log10(abs(fftshift(1/sqrt(len)*fft(x_base_t)))))
% xlabel('DTFT frequency f in MHz')
% ylabel('|X^{base}(f)| in dB')
% title('Frequency Domain Plot for transmitted signal, x_{base}(t)')
% xline(11, 'red');
% xline(-11, 'red');
% xline(15, 'black');
% xline(-15, 'black');
% yline(-20, 'red');
% yline(-40, 'black');
% axis tight

%% Received Signal - y_base(t)

y_base_t = receivedsignal;
y_base_t = y_base_t(1:30000);
len = length(y_base_t);
t_microseconds=[0:len-1]/200e6*1e6;

figure(2)
LargeFigure(gcf, 0.15); % Make figure large
clf

% Plot time domain signal
subplot(2,1,1);
hold on;
plot(t_microseconds, real(y_base_t))
plot(t_microseconds, imag(y_base_t))
hold off;
title('Time Domain Plot for received signal, y_{base}(t)')
xlabel('t in microseconds')
ylabel('y^{base}(t)')
legend("real", "imag");
axis tight

% Plot frequency domain signal
subplot(2,1,2);
plot([-len/2+1:len/2]*200/len, 20*log10(abs(fftshift(1/sqrt(len)*fft(y_base_t)))))
xlabel('DTFT frequency f in MHz')
ylabel('|Y^{base}(f)| in dB')
title('Frequency Domain Plot for received signal, y_{base}(t)')
xline(11, 'red');
xline(-11, 'red');
xline(15, 'black');
xline(-15, 'black');
yline(-20, 'red');
yline(-40, 'black');
axis tight


% %% Pulse - p^b(t)
% 
% len = length(pulse);
% t_microseconds=[0:len-1]/200e6*1e6;
% 
% figure(3)
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% 
% % Plot time domain signal
% subplot(2,1,1);
% hold on;
% plot(t_microseconds, real(pulse))
% plot(t_microseconds, imag(pulse))
% hold off;
% title('Time Domain Plot for Pulse')
% xlabel('t in microseconds')
% ylabel('p^{b}(t)')
% legend("real", "imag");
% axis tight
% 
% % Plot frequency domain signal
% subplot(2,1,2);
% plot([-len/2+1:len/2]*200/len, 20*log10(abs(fftshift(1/sqrt(len)*fft(pulse)))))
% xlabel('DTFT frequency f in MHz')
% ylabel('|P^{b}(f)| in dB')
% title('Frequency Domain Plot for Pulse')
% axis tight

%% Sampler Output - z_{k}
len = length(zk);
samples=[0:len-1];

figure(4)
% LargeFigure(gcf, 0.15); % Make figure large
clf

% Plot time domain signal
hold on;
plot(samples, real(zk))
plot(samples, imag(zk))
hold off;
title('Sampler Output (z_{k})')
xlabel('Samples')
ylabel('z_{k}')
legend("real", "imag");
axis tight

%% y_base(t) Post-Time Recovery

len = length(zt_inphase_timing_plot);
t_microseconds=[0:len-1]/200e6*1e6;

figure(5)
% LargeFigure(gcf, 0.15); % Make figure large
clf

% Plot time domain signal
hold on;
plot(t_microseconds, real(zt_inphase_timing_plot))
plot(t_microseconds, imag(zt_inphase_timing_plot))
hold off;
title('y^{base}(t) after Time Recovery')
xlabel('t in microseconds')
ylabel('y^{base}(t)')
legend("real", "imag");
axis tight

% %% Equalizer sample output - v_k
% 
% vk = [zk1; zk2; zk3; zk4];
% vk = vk(1:1500);
% len = length(vk);
% samples=[0:len-1];
% 
% figure(6)
% % LargeFigure(gcf, 0.15); % Make figure large
% clf
% 
% % Plot time domain signal
% hold on;
% plot(samples, real(vk))
% plot(samples, imag(vk))
% hold off;
% title('Equalizer Sample Output (v_{k})')
% xlabel('Samples')
% ylabel('v_{k}')
% legend("real", "imag");
% axis tight
% 
% figure(8)
% hold on;
% scatter(real(zk1), imag(zk1));
% scatter(real(zk2), imag(zk2));
% scatter(real(zk3), imag(zk3));
% scatter(real(zk4(1:240)), imag(zk4(1:240)));
% scatter([-1, 1], [0,0]);
