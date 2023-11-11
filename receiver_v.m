clear; close all; clc;
load freq_sync;
load t_sync.mat;
load frame_sync.mat;
load transmitsignal_test_cases.mat;
load receivedsignal_test_cases.mat;
sync_len = 1; % microseconds
fs = 200; %MHz
upsampling_rate = 20;
pulse = rcosdesign(0.5, 7, upsampling_rate, 'sqrt');

preamble = horzcat(freq_sync, t_sync, frame_sync);
impulse = horzcat(zeros(1, 100), 1, zeros(1, 100));
sine = horzcat(zeros(1, 100), sin(2*pi*5e6.*linspace(0, 50e-6, 500)), zeros(1, 100));
square = horzcat(zeros(1, 100), ones(1, 100), zeros(1, 100));

%% Timing Synchronization

%Demodulate
wt = fliplr(pulse); %probably unnecessary
zt = conv(wt, receivedsignal);

sig_to_find = t_sync;  %PUT SIGNAL TO FIND HERE

xt_conj = conj(sig_to_find); 

tau = length(zt);
sums = zeros(1, tau);
p = length(xt_conj);

for i = 1:tau-p*upsampling_rate-1
    sums(i) = 0;
    for k = 1:p
        sums(i) = sums(i) + xt_conj(k)*zt(k*upsampling_rate+i);
    end
end

timing_offset = find(sums == max(sums));
zt = zt(timing_offset:end);
phase = angle(max(sums));
zt_inphase_timing = zt .* exp(-j*phase); % Remove phase offset

%% Frame Synchronization

sig_to_find = frame_sync;  %PUT SIGNAL TO FIND HERE

xt_conj = conj(sig_to_find);

tau = length(zt);
sums = zeros(1, tau);
p = length(xt_conj);

for i = 1:tau-p*upsampling_rate-1
    sums(i) = 0;
    for k = 1:p
        sums(i) = sums(i) + xt_conj(k)*zt(k*upsampling_rate+i);
    end
end

offset = find(sums == max(sums));
phase = angle(max(sums));

zt_inphase_frame = zt .* exp(-j*phase); % Remove phase offset

start_ind = offset + length(frame_sync);
zt_inphase_frame = zt_inphase_frame(start_ind:end);

%% Equalization
received_preamble = zt_inphase_frame(start_ind:start_ind+length(preamble)-1);
channel_effect = received_preamble./preamble;
h0=mean(channel_effect);
% h0=0.2993;
zt_inphase_frame = zt_inphase_frame./h0;

%% Plotting Results
xt = zt_inphase_frame;
len = length(xt);
% Plot time domain signal
% subplot(2,1,1);
hold on;
plot([0:len-1]/200e6*1e6, real(xt))
plot([0:len-1]/200e6*1e6, imag(xt))
hold off;
xlabel('t in microseconds')
ylabel('x(t)')
legend("real", "imag");
title('xbase(t)')
axis tight
