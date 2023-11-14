%% Setup
clear; clc; close all;

sync_len = 5; % microseconds
fs = 200; %MHz
upsampling_rate = 20;

% freq_sync = ones(1, fs*sync_len/upsampling_rate); %2 us of frame sync for the receiver to get PLL lock
% t_sync = floor(2.*rand(1, fs*sync_len/upsampling_rate)).* 2 -1; %Pseudo-random BPSK
% frame_sync = floor(2.*rand(1, fs*sync_len/upsampling_rate)).* 2 -1; %Pseudo-random BPSK
% save("frame_sync", "frame_sync")
% save("t_sync", "t_sync")
% save("freq_sync", "freq_sync")

load freq_sync;
load t_sync.mat;
load frame_sync.mat;

preamble = horzcat(freq_sync, t_sync, frame_sync);
% impulse = horzcat(zeros(1, 100), 1, zeros(1, 100));
% sine = horzcat(zeros(1, 100), sin(2*pi*5e6.*linspace(0, 50e-6, 500)), zeros(1, 100));
% square = horzcat(zeros(1, 100), ones(1, 100), zeros(1, 100));
% test_sig = horzcat(preamble, impulse, zeros(1,40), sine, zeros(1, 40), square);
% test_sig = upsample(test_sig, upsampling_rate);

data = imread('images/shannon1440.bmp');
dims = size(data);
bits = reshape(data, 1, dims(1)*dims(2));
bits = bits .* 2 -1; %BPSK
data1 = bits(1:400);
data2 = bits(401:800);
data3 = bits(801:1200);
data4 = bits(1201:min(length(bits), 1600));
bits = horzcat(preamble, data1, frame_sync, data2, frame_sync, data3, frame_sync, data4);
bits_up = upsample(bits, upsampling_rate);
pulse = rcosdesign(0.9, 30, upsampling_rate, 'sqrt');
save('bits', 'bits');

transmitsignal = conv(bits_up, pulse);
save("transmitsignal.mat", "transmitsignal");
time = length(transmitsignal)/200; %In microseconds

