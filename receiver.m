clear; close all; clc;
load freq_sync;
load qam_bits;
load t_sync.mat;
load frame_sync.mat;
load receivedsignal.mat;
load transmitsignal.mat;
load bits.mat

sigman = 0.2;

transmitsignalwithdelay = [zeros(1, 2147), transmitsignal];
receivedsignal = exp(j*pi/6) * transmitsignalwithdelay + sigman/sqrt(2) * (randn(size(transmitsignalwithdelay))+j*randn(size(transmitsignalwithdelay)));
receivedsignal=receivedsignal';

sync_len = 1; % microseconds
fs = 200; %MHz
M = 16;
b=log2(M);
upsampling_rate = 12;
pulse = rcosdesign(0.9, 30, upsampling_rate, 'sqrt');

%16-QAM constellation
d = sqrt(2)/3;
options = [1.5+1.5j, 0.5+1.5j, -1.5+1.5j, -0.5+1.5j
           1.5+0.5j, 0.5+0.5j, -1.5+0.5j, -0.5+0.5j
           1.5-1.5j, 0.5-1.5j, -1.5-1.5j, -0.5-1.5j
           1.5-0.5j, 0.5-0.5j, -1.5-0.5j, -0.5-0.5j].*d;

preamble = horzcat(freq_sync, t_sync, frame_sync);
impulse = horzcat(zeros(1, 100), 1, zeros(1, 100));
sine = horzcat(zeros(1, 100), sin(2*pi*5e6.*linspace(0, 50e-6, 500)), zeros(1, 100));
square = horzcat(zeros(1, 100), ones(1, 100), zeros(1, 100));

%% Filter
wt = fliplr(pulse); %probably unnecessary
zt = conv(receivedsignal, wt);

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
for i = 1:num_frames
   phases(i) = angle(sums(locs(i)));
   starts(i) = locs(i) + upsampling_rate*((length(frame_sync)/b)+.55);
end

%% Equalization
received_preamble = zt_inphase_timing(1:upsampling_rate:length(preamble)*upsampling_rate);
channel_effect = abs(received_preamble)'./abs(preamble);
h0=mean(abs(channel_effect));
zt_inphase_timing = zt_inphase_timing./h0;

%% Sampling
% zk = zt_inphase_frame(1:upsampling_rate:length(zt_inphase_frame));
% num_frames = 10;
% 
% frames = zeros( num_frames);
% for i = 1:num_frames
%     i
%     frames(i) = zt_inphase_timing(starts(i): starts(i+1)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(i));
% end
% 
% zks = zeros(1, num_frames);
% for i = 1:num_frames
%     zks(i) = zt_inphase_timing(starts(i): starts(i+1)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(i));
% end
% 
% z_demod = qamdemod(zks, M);
% 


frame1 = zt_inphase_timing(starts(1): starts(2)-upsampling_rate*(length(frame_sync)/b)).* exp(-j*phases(1));
frame2 = zt_inphase_timing(starts(2): starts(3)-upsampling_rate*(length(frame_sync)/b)).* exp(-j*phases(2));
frame3 = zt_inphase_timing(starts(3): starts(4)-upsampling_rate*(length(frame_sync)/b)).* exp(-j*phases(3));
frame4 = zt_inphase_timing(starts(4): starts(5)-upsampling_rate*(length(frame_sync)/b)).* exp(-j*phases(4));
frame5 = zt_inphase_timing(starts(5): starts(6)-upsampling_rate*(length(frame_sync)/b)).* exp(-j*phases(5));
frame6 = zt_inphase_timing(starts(6): starts(7)-upsampling_rate*(length(frame_sync)/b)).* exp(-j*phases(6));
frame7 = zt_inphase_timing(starts(7): starts(8)-upsampling_rate*(length(frame_sync)/b)).* exp(-j*phases(7));
frame8 = zt_inphase_timing(starts(8): starts(9)-upsampling_rate*(length(frame_sync)/b)).* exp(-j*phases(8));
frame9 = zt_inphase_timing(starts(9): starts(10)-upsampling_rate*(length(frame_sync)/b)).* exp(-j*phases(9));
frame10 = zt_inphase_timing(starts(10): end-1-upsampling_rate*(length(frame_sync)/b)).* exp(-j*phases(10));

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

zk = [zk1; zk2; zk3; zk4; zk5; zk6; zk7; zk8; zk9; zk10; ];

zk_demod = qamdemod(zk, M);
len = length(zk_demod);

% Find the closest symbol in the constellation for each sample
rx_symbols = zeros(1,len*b);

for i = 1:len
    [~, index] = min(abs(zk_demod(i) - options(:)));
    rx_symbols(((i-1)*b)+1: i*b) = de2bi(index-1,b);
end



%% BER
message = imread("images/shannon1440.bmp");
message_vec = reshape(message, 1, []);

bits = message_vec;

BER = mean(rx_symbols ~= bits)
disp(['BER is ', num2str(BER)])

% % Guessing
% rx_ipts1 = sign(real(zk1));
% rx_ipts2 = sign(real(zk2));
% rx_ipts3 = sign(real(zk3));
% rx_ipts4 = sign(real(zk4));
% 
% % Error checking and BER calculation
% data = imread('images/shannon1440.bmp');
% dims = size(data);
% bits = reshape(data, 1, dims(1)*dims(2));
% bits = bits .* 2 -1; %BPSK
% data1 = bits(1:400);
% data2 = bits(401:800);
% data3 = bits(801:1200);
% data4 = bits(1201:min(length(bits), 1600));
% 
% BER1 = length(find(rx_ipts1(1:400) ~= data1'))/length(data1);
% BER2 = length(find(rx_ipts2(1:400) ~= data2'))/length(data2);
% BER3 = length(find(rx_ipts3(1:400) ~= data3'))/length(data3);
% BER4 = length(find(rx_ipts4(1:240) ~= data4'))/length(data4);
% 
% BER = ((BER1+BER2+BER3)*400+BER4*240)/1440

%% Plotting

figure; imshow(message); title("original")

% received_image = [rx_ipts1(1:400); rx_ipts2(1:400); rx_ipts3(1:400); rx_ipts4(1:240)];
% received_image = reshape(received_image, 45, 32);
figure; imshow(zk_demod); 


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
% 
% %% Received Signal - y_base(t)
% 
% y_base_t = receivedsignal;
% y_base_t = y_base_t(1:30000);
% len = length(y_base_t);
% t_microseconds=[0:len-1]/200e6*1e6;
% 
% figure(2)
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% 
% % Plot time domain signal
% subplot(2,1,1);
% hold on;
% plot(t_microseconds, real(y_base_t))
% plot(t_microseconds, imag(y_base_t))
% hold off;
% title('Time Domain Plot for received signal, y_{base}(t)')
% xlabel('t in microseconds')
% ylabel('y^{base}(t)')
% legend("real", "imag");
% axis tight
% 
% % Plot frequency domain signal
% subplot(2,1,2);
% plot([-len/2+1:len/2]*200/len, 20*log10(abs(fftshift(1/sqrt(len)*fft(y_base_t)))))
% xlabel('DTFT frequency f in MHz')
% ylabel('|Y^{base}(f)| in dB')
% title('Frequency Domain Plot for received signal, y_{base}(t)')
% xline(11, 'red');
% xline(-11, 'red');
% xline(15, 'black');
% xline(-15, 'black');
% yline(-20, 'red');
% yline(-40, 'black');
% axis tight
% 
% 
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
% 
% %% Sampler Output - z_{k}
% frame1 = zt_inphase_timing_plot(starts(1): starts(2)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(1));
% frame2 = zt_inphase_timing_plot(starts(2): starts(3)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(2));
% frame3 = zt_inphase_timing_plot(starts(3): starts(4)-upsampling_rate*(length(frame_sync))).* exp(-j*phases(3));
% frame4 = zt_inphase_timing_plot(starts(4): end-1-upsampling_rate*(length(frame_sync))).* exp(-j*phases(4));
% 
% z1 = frame1(1:upsampling_rate:length(frame1));
% z2 = frame2(1:upsampling_rate:length(frame2));
% z3 = frame3(1:upsampling_rate:length(frame3));
% z4 = frame4(1:upsampling_rate:length(frame4));
% 
% zk = [z1; z2; z3; z4];
% len = length(zk);
% samples=[0:len-1];
% 
% figure(4)
% % LargeFigure(gcf, 0.15); % Make figure large
% clf
% 
% % Plot time domain signal
% hold on;
% plot(samples, real(zk))
% plot(samples, imag(zk))
% hold off;
% title('Sampler Output (z_{k})')
% xlabel('Samples')
% ylabel('z_{k}')
% legend("real", "imag");
% axis tight
% 
% %% y_base(t) Post-Time Recovery
% 
% len = length(zt_inphase_timing_plot);
% t_microseconds=[0:len-1]/200e6*1e6;
% 
% figure(5)
% % LargeFigure(gcf, 0.15); % Make figure large
% clf
% 
% % Plot time domain signal
% hold on;
% plot(t_microseconds, real(zt_inphase_timing_plot))
% plot(t_microseconds, imag(zt_inphase_timing_plot))
% hold off;
% title('y^{base}(t) after Time Recovery')
% xlabel('t in microseconds')
% ylabel('y^{base}(t)')
% legend("real", "imag");
% axis tight
% 
% %% Equalizer sample output - v_k
% 
% vk = [zk1; zk2; zk7; zk4];
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
