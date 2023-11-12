clear; close all; clc;
load freq_sync;
load t_sync.mat;
load frame_sync.mat;
load receivedsignal.mat;
load transmitsignal.mat;
load bits.mat
sync_len = 1; % microseconds
fs = 200; %MHz
upsampling_rate = 20;
pulse = rcosdesign(0.9, 30, upsampling_rate, 'sqrt');

preamble = horzcat(freq_sync, t_sync, frame_sync);
impulse = horzcat(zeros(1, 100), 1, zeros(1, 100));
sine = horzcat(zeros(1, 100), sin(2*pi*5e6.*linspace(0, 50e-6, 500)), zeros(1, 100));
square = horzcat(zeros(1, 100), ones(1, 100), zeros(1, 100));

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

%% Frame Synchronization

sig_to_find = frame_sync;  %PUT SIGNAL TO FIND HERE
num_frames = 4;
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

[pks,locs] = findpeaks(abs(sums), 'MinPeakProminence', 10);
phases = zeros(1, num_frames);
starts = zeros(1, num_frames);
for i = 1:4
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
frame4 = zt_inphase_timing(starts(4): end-1-upsampling_rate*(length(frame_sync))).* exp(-j*phases(4));

zk1 = frame1(1:upsampling_rate:length(frame1));
zk2 = frame2(1:upsampling_rate:length(frame2));
zk3 = frame3(1:upsampling_rate:length(frame3));
zk4 = frame4(1:upsampling_rate:length(frame4));

% Guessing
rx_ipts1 = sign(real(zk1));
rx_ipts2 = sign(real(zk2));
rx_ipts3 = sign(real(zk3));
rx_ipts4 = sign(real(zk4));

%Error checking and BER calculation
data = imread('images/shannon1440.bmp');
dims = size(data);
bits = reshape(data, 1, dims(1)*dims(2));
bits = bits .* 2 -1; %BPSK
data1 = bits(1:400);
data2 = bits(401:800);
data3 = bits(801:1200);
data4 = bits(1201:min(length(bits), 1600));

BER1 = length(find(rx_ipts1(1:400) ~= data1'))/length(data1)
BER2 = length(find(rx_ipts2(1:400) ~= data2'))/length(data2)
BER3 = length(find(rx_ipts3(1:400) ~= data3'))/length(data3)
BER4 = length(find(rx_ipts4(1:240) ~= data4'))/length(data4)

%% Plotting

figure; imshow(data); title("original")

received_image = [rx_ipts1(1:400); rx_ipts2(1:400); rx_ipts3(1:400); rx_ipts4(1:240)];
received_image = reshape(received_image, 45, 32);
figure; imshow(received_image); 


% xt = upsample(data4, upsampling_rate);
% len = length(xt);
% 
% figure(1)
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% 
% % Plot transmit x(t)
% subplot(2,1,1);
% hold on;
% plot([0:len-1], real(xt))
% hold off;
% xlabel('t in microseconds')
% ylabel('f(t)')
% title('xbase(t)')
% axis tight
% 
% xt = zk4;
% len = length(xt);
% 
% % Plot received signal of choice
% subplot(2,1,2);
% hold on
% plot([0:len-1], real(xt))
% plot([0:len-1], imag(xt))
% hold off;
% title('z(t) frame in phase')
% xlabel('t in microseconds')
% ylabel('x(t)')
% legend("real", "imag");
% axis tight
