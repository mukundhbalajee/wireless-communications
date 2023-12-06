clear; close all; clc;
load receivedsignal2.mat;
load transmitsignal.mat;
load qam_pre.mat
load qam_fsync.mat

MMSE = true;

preamble = qam_preamble;
frame_sync = qam_frame_sync;

sync_len = 12;
% microseconds
fs = 200; %MHz
upsampling_rate = 12;
pulse = rcosdesign(0.3, 30, upsampling_rate, 'sqrt');

%% Demodulate
wt = fliplr(pulse); %probably unnecessary
zt = conv(wt, receivedsignal);

%% Timing Synchronization

num_frames = 20;

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

[pks,inds] = findpeaks(abs(sums), 'MinPeakProminence', 5);
maxes = maxk(pks, num_frames);

timing_offset = inds(1);
zt_timing = zt(timing_offset:end);
phase = angle(sums(inds(1)));
zt_inphase_timing = zt_timing;
zt_inphase_timing_plot = zt_inphase_timing;

%% Frame Synchronization

sig_to_find = frame_sync;  %PUT SIGNAL TO FIND HERE
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
starts_pilot = locs;
phases = angle(sums(locs));
for i = 1:num_frames-1
    if(phases(i+1) < phases(i))
        phases(i+1:end) = phases(i+1:end) + 2*pi;
    end
end
starts = locs + upsampling_rate*(length(frame_sync)+1);

%% Equalization

if(~MMSE)
    disp("One Tap EQ");
    received_preamble = zt_inphase_timing(1:upsampling_rate:length(preamble)*upsampling_rate);
    channel_effect = abs(received_preamble)'./abs(preamble);
    h0=mean(abs(channel_effect));
    zt_inphase_timing = zt_inphase_timing./h0;
end

%% Sampling
frame1 = zt_inphase_timing(starts(1): starts(2)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(1)+phases(2)));
frame2 = zt_inphase_timing(starts(2): starts(3)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(2)+phases(3)));
frame3 = zt_inphase_timing(starts(3): starts(4)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(3)+phases(4)));
frame4 = zt_inphase_timing(starts(4): starts(5)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(4)+phases(5)));
frame5 = zt_inphase_timing(starts(5): starts(6)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(5)+phases(6)));
frame6 = zt_inphase_timing(starts(6): starts(7)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(6)+phases(7)));
frame7 = zt_inphase_timing(starts(7): starts(8)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(7)+phases(8)));
frame8 = zt_inphase_timing(starts(8): starts(9)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(8)+phases(9)));
frame9 = zt_inphase_timing(starts(9): starts(10)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(9)+phases(10)));
frame10 = zt_inphase_timing(starts(10): starts(11)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(10)+phases(11)));
frame11 = zt_inphase_timing(starts(11): starts(12)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(11)+phases(12)));
frame12 = zt_inphase_timing(starts(12): starts(13)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(12)+phases(13)));
frame13 = zt_inphase_timing(starts(13): starts(14)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(13)+phases(14)));
frame14 = zt_inphase_timing(starts(14): starts(15)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(14)+phases(15)));
frame15 = zt_inphase_timing(starts(15): starts(16)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(15)+phases(16)));
frame16 = zt_inphase_timing(starts(16): starts(17)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(16)+phases(17)));
frame17 = zt_inphase_timing(starts(17): starts(18)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(17)+phases(18)));
frame18 = zt_inphase_timing(starts(18): starts(19)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(18)+phases(19)));
frame19 = zt_inphase_timing(starts(19): starts(20)-upsampling_rate*(length(frame_sync))).* exp(-j*0.5*(phases(19)+phases(20)));
frame20 = zt_inphase_timing(starts(20): end).* exp(-j*phases(20)+0.5*(phases(20)-phases(19)));

len = upsampling_rate*floor(length(frame1)/upsampling_rate);
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
zk11 = frame11(1:upsampling_rate:len);
zk12 = frame12(1:upsampling_rate:len);
zk13 = frame13(1:upsampling_rate:len);
zk14 = frame14(1:upsampling_rate:len);
zk15 = frame15(1:upsampling_rate:len);
zk16 = frame16(1:upsampling_rate:len);
zk17 = frame17(1:upsampling_rate:len);
zk18 = frame18(1:upsampling_rate:len);
zk19 = frame19(1:upsampling_rate:len);
zk20 = frame20(1:upsampling_rate:len);

zk = [zk1; zk2; zk3; zk4; zk5; zk6; zk7; zk8; zk9; zk10; zk11; zk12; zk13; zk14; zk15; zk16; zk17; zk18; zk19; zk20];
zk = transpose(zk);

%% Equalization - MMSE-LE
if(MMSE)
    disp("MMSE-LE EQ");
    filter_length = 6;
    mu = .2;
    eq_weights = zeros(filter_length, num_frames);
    
    %Train filter based on each received chunk
    for frame = 1:num_frames
        received_pilot_frame = zt_inphase_timing(starts_pilot(frame) + upsampling_rate:(length(frame_sync)+1)*upsampling_rate+starts_pilot(frame)).* exp(-j*phases(frame));
        received_pilot = received_pilot_frame(1:upsampling_rate:length(received_pilot_frame));
        trained_filter = zeros(1, filter_length);
        trained_filter = train_filter(trained_filter, 100, received_pilot, frame_sync, mu, filter_length);
        eq_weights(:,frame) = trained_filter.';
    end
    
    zk_eq = zeros(size(zk));
    frame_size = (length(zk) / num_frames);
    start_ind = 1;
    for frame = 1:num_frames
        end_ind = start_ind+frame_size-1;
        if frame == num_frames
            end_ind = length(zk);
        end
        frame_size = end_ind-start_ind+1;
        received_frame = zk(start_ind:end_ind);
        eq_val = conv(received_frame, eq_weights(:, frame));
        zk_eq(start_ind:end_ind) = eq_val(1:frame_size);
        start_ind = end_ind+1;
    end
    zk_eq = transpose(zk_eq);
    zk = zk_eq;
end

%% Guessing
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

%% QAM to Bits
rx_bits = dec2bin(received_pts, 4);
rx_bits = split(char(strjoin(string(rx_bits),'')),'');
rx_bits = rx_bits(2:end-1);
rx_bits = str2num(cell2mat(rx_bits));
rx_bits = rx_bits';

%% Error checking and BER calculation
image = imread('images/shannon20520.bmp');
dims = size(image);

rx_bits = rx_bits(1:dims(1)*dims(2));
bits = double(reshape(image, 1, dims(1)*dims(2)));

start = 1;
stop = length(bits)/20*20;
BER = length(find(rx_bits(start:stop) ~= bits(start:stop)))/(stop-start)

%% Plotting

received_image = [rx_bits];
received_image = reshape(received_image, dims(1), dims(2));

figure(10); 

subplot(2,1,1);
hold on;
imshow(image);  title("original")
hold off;

received_image = [rx_bits];
received_image = reshape(received_image, dims(1), dims(2));

subplot(2,1,2); 
hold on;
imshow(received_image); title("received")
hold off;


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

%% Received Signal - y_base(t)

y_base_t = receivedsignal;;
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


%% Pulse - p^b(t)

len = length(pulse);
t_microseconds=[0:len-1]/200e6*1e6;

figure(3)
LargeFigure(gcf, 0.15); % Make figure large
clf

% Plot time domain signal
subplot(2,1,1);
hold on;
plot(t_microseconds, real(pulse))
plot(t_microseconds, imag(pulse))
hold off;
title('Time Domain Plot for Pulse')
xlabel('t in microseconds')
ylabel('p^{b}(t)')
legend("real", "imag");
axis tight

% Plot frequency domain signal
subplot(2,1,2);
plot([-len/2+1:len/2]*200/len, 20*log10(abs(fftshift(1/sqrt(len)*fft(pulse)))))
xlabel('DTFT frequency f in MHz')
ylabel('|P^{b}(f)| in dB')
title('Frequency Domain Plot for Pulse')
axis tight

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

%% Equalizer sample output - v_k

vk = zk;
len = length(vk);
samples=[0:len-1];

figure(6)
% LargeFigure(gcf, 0.15); % Make figure large
clf

% Plot time domain signal
hold on;
plot(samples, real(vk))
plot(samples, imag(vk))
hold off;
title('Equalizer Sample Output (v_{k})')
xlabel('Samples')
ylabel('v_{k}')
legend("real", "imag");
axis tight
%% Constellation Plot

figure(8)
hold on;
scatter(real(zk), imag(zk));
scatter(reshape(real(options), 1, 16), reshape(imag(options), 1, 16), 'LineWidth', 2);
xlabel('I'); ylabel('Q')
hold off;

%% LMS Algorithm

function [eq_weights] = train_filter(eq_weights, train_iter, received_signal, expected_signal, mu, filter_length)
    for j = 1:train_iter
        for i = filter_length:length(expected_signal)
            received_segment = received_signal(i:-1:i-filter_length+1);
            predicted_output = eq_weights * received_segment;
            error = predicted_output - expected_signal(i);
            new_w = mu * error * conj(received_segment).';
            eq_weights = eq_weights - reshape(new_w, 1, filter_length);
        end
    end
end
