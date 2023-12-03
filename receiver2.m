clear; close all; clc;
load receivedsignal.mat;
load transmitsignal.mat;
load qam_pre.mat
load qam_fsync.mat

preamble = qam_preamble;
frame_sync = qam_frame_sync;

% receivedsignal = transmitsignal.';

sync_len = 24; % microseconds
fs = 200; %MHz
upsampling_rate = 12;
M = 16;
b = log2(M);
pulse = rcosdesign(0.3, 30, upsampling_rate, 'sqrt');

%% Matched Filter
wt = fliplr(pulse); %probably unnecessary
zt = conv(wt, receivedsignal);

%% Timing Synchronization

sig_to_find = preamble;

xt_conj = conj(sig_to_find); 

tau = length(zt);
sums = zeros(1, tau);
p = length(xt_conj);

for i = 1:tau-p*upsampling_rate
    sums(i) = 0;
    for k = 1:p
        sums(i) = sums(i) + xt_conj(k)*zt((k-1)*upsampling_rate+i);
    end
end

timing_offset = find(abs(sums) == max(abs(sums)));
zt_timing = zt(timing_offset:end);
phase = angle(max(sums));
zt_inphase_timing = zt_timing .* exp(-j*phase); % Remove phase offset
zt_inphase_timing_plot = zt_inphase_timing;

%% Frame Synchronization

sig_to_find = frame_sync;
num_frames = 10;
xt_conj = conj(sig_to_find);

tau = length(zt_inphase_timing);
sums = zeros(1, tau);
p = length(xt_conj);

for i = 1:tau-p*upsampling_rate-1
    sums(i) = 0;
    for k = 1:p
        sums(i) = sums(i) + xt_conj(k)*zt_inphase_timing((k-1)*upsampling_rate+i);
    end
end

[pks,inds] = findpeaks(abs(sums), 'MinPeakProminence', 1);
maxes = maxk(pks, num_frames);

locs = zeros(1, num_frames);

for i = 1:num_frames
    locs(i) = find(abs(sums) == maxes(i));
end

locs = sort(locs);
starts_pilot = locs;
phases = angle(sums(locs));
starts = locs + upsampling_rate*(length(frame_sync));

%% Sampling
frame1 = zt_inphase_timing(starts(1): starts(2)-upsampling_rate*(length(frame_sync)));
frame2 = zt_inphase_timing(starts(2): starts(3)-upsampling_rate*(length(frame_sync)));
frame3 = zt_inphase_timing(starts(3): starts(4)-upsampling_rate*(length(frame_sync)));
frame4 = zt_inphase_timing(starts(4): starts(5)-upsampling_rate*(length(frame_sync)));
frame5 = zt_inphase_timing(starts(5): starts(6)-upsampling_rate*(length(frame_sync)));
frame6 = zt_inphase_timing(starts(6): starts(7)-upsampling_rate*(length(frame_sync)));
frame7 = zt_inphase_timing(starts(7): starts(8)-upsampling_rate*(length(frame_sync)));
frame8 = zt_inphase_timing(starts(8): starts(9)-upsampling_rate*(length(frame_sync)));
frame9 = zt_inphase_timing(starts(9): starts(10)-upsampling_rate*(length(frame_sync)));
frame10 = zt_inphase_timing(starts(10): end);

len = 513*upsampling_rate;
zk1 = frame1(1:upsampling_rate:len);
zk2 = frame2(1:upsampling_rate:len);
zk3 = frame3(1:upsampling_rate:len);
zk4 = frame4(1:upsampling_rate:len);
zk5 = frame5(1:upsampling_rate:len);
zk6 = frame6(1:upsampling_rate:len);
zk7 = frame7(1:upsampling_rate:len);
zk8 = frame8(1:upsampling_rate:len);
zk9 = frame9(1:upsampling_rate:len);
zk10 = frame10(1:upsampling_rate:len);

zk = [zk1; zk2; zk3; zk4; zk5; zk6; zk7; zk8; zk9; zk10];
zk = transpose(zk);

%% Equalization - MMSE-LE

filter_length = 21;
mu = .3;
eq_weights = zeros(filter_length, num_frames);

%Train filter based on each received chunk                     
for frame = 1:num_frames
    received_pilot_frame = zt_inphase_timing(starts_pilot(frame):(length(frame_sync)-1)*upsampling_rate+starts_pilot(frame));
    received_pilot = received_pilot_frame(1:upsampling_rate:length(received_pilot_frame));
    trained_filter = eq_weights(:, frame);
    trained_filter = train_filter(trained_filter, 50, received_pilot, frame_sync, mu, filter_length);
    eq_weights(:,frame) = trained_filter;
end

zk_eq = zeros(size(zk));
frame_size = (length(zk) / num_frames);
start_ind = 1;
for frame = 1:num_frames
    end_ind = start_ind+frame_size-1;
    if frame == num_frames
        end_ind = length(zk);
    end
    received_frame = zk(start_ind:end_ind);
    eq_val = conv(received_frame, eq_weights(:, frame));
    zk_eq(start_ind:end_ind) = eq_val(1:frame_size);
    start_ind = end_ind+1;
end

%% Guessing
d = sqrt(2)/3;
options = [1.5+1.5j, 0.5+1.5j, -1.5+1.5j, -0.5+1.5j
           1.5+0.5j, 0.5+0.5j, -1.5+0.5j, -0.5+0.5j
           1.5-1.5j, 0.5-1.5j, -1.5-1.5j, -0.5-1.5j
           1.5-0.5j, 0.5-0.5j, -1.5-0.5j, -0.5-0.5j].*d;

received_pts = zeros(1, length(zk_eq)); %QAM

min_pt = -1;
min_distance = 10000;
for i = 1:length(zk_eq)
    min_pt = -1;
    min_distance = 10000;
    for j = 1:M
        refpoint = [real(options(j)), imag(options(j))];
        point = [real(zk_eq(i)), imag(zk_eq(i))];
        distance = sqrt((refpoint(1)-point(1))^2 + (refpoint(2)-point(2))^2);
        if(distance < min_distance)
            min_distance = distance;
            min_pt = j;
        end
    end
    received_pts(i) = min_pt-1;
end

%% QAM to Bits
rx_bits = dec2bin(received_pts, 4);
rx_bits = split(char(strjoin(string(rx_bits),'')),'');
rx_bits = rx_bits(2:end-1);
rx_bits = str2num(cell2mat(rx_bits));
rx_bits = rx_bits';

%% Error checking and BER calculation
image = imread('images/shannon20520.bmp');
dims = size(image);

bits = string(double(reshape(image, 4, dims(1)*dims(2)/4)));

qam_bits = strings(1, length(bits));
for i = 1:length(bits)
    qam_bits(i) = strjoin(bits(:,i),'');
end

qam_bits = bin2dec((qam_bits))+1;
qam_points = options(qam_bits);

bits = double(reshape(image, 1, dims(1)*dims(2)));
n=8;
start = 1;%n*4*length(zk1)+1;
stop = length(bits); %start + 400; %(n+1)*4*length(zk1);
BER = length(find(rx_bits(start:stop) ~= bits(start:stop)))/(stop-start)*100
N = 1500;
real_diffs = (real((qam_points(1:N))) - real(zk_eq(1:N)));
imag_diffs = (imag((qam_points(1:N))) - imag(zk_eq(1:N)));
mean(real_diffs)
mean(imag_diffs)

%% Plotting

% figure(1); 
% 
% subplot(2,1,1);
% hold on;
% imshow(image); title("original")
% hold off;
% 
% received_image = [rx_bits];
% received_image = reshape(received_image, dims(1), dims(2));
% 
% subplot(2,1,2); 
% hold on;
% imshow(received_image); title("received")
% hold off;

%% Transmitted signal - x_base(t)

x_base_t = transmitsignal;
len = length(x_base_t);
t_microseconds=[0:len-1]/200e6*1e6;

figure(2);
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

vk = zk_eq;
len = length(vk);
samples=[0:len-1];

subplot(2,1,2);
hold on;
plot(samples, real(vk))
plot(samples, imag(vk))
hold off;
title('Equalizer Sample Output (v_{k})')
xlabel('Samples')
ylabel('v_{k}')
legend("real", "imag");
axis tight

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
% figure(3);
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
% figure(4);
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
% len = length(zk);
% samples=[0:len-1];
% 
% figure(5);
% LargeFigure(gcf, 0.15); % Make figure large
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
% figure(6);
% LargeFigure(gcf, 0.15); % Make figure large
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

% %% Equalizer sample output - v_k
% 
% vk = zk_eq;
% len = length(vk);
% samples=[0:len-1];
% 
% figure(7);
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% 
% subplot(2,1,1);
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
% subplot(2,1,2);
% hold on;
% scatter(real(zk_eq), imag(zk_eq));
% scatter([-1, 1], [0,0]);
% legend("real", "imag");
% axis tight

%% LMS Algorithm
function [eq_weights] = train_filter(eq_weights, train_iter, received_signal, expected_signal, mu, filter_length)
    for j = 1:train_iter
        for i = filter_length:length(expected_signal)
            received_segment = received_signal(i:-1:i-filter_length+1);
            predicted_output = eq_weights .* received_segment;
            error = predicted_output - expected_signal(i);
            new_w = mu * error .* conj(received_segment);
            eq_weights = eq_weights - new_w;
        end
    end
end
