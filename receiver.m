load freq_sync;
load t_sync.mat;
load frame_sync.mat;
load transmitsignal.mat;
load receivedsignal.mat;
sync_len = 1; % microseconds
fs = 200; %MHz
upsampling_rate = 20;
pulse = rcosdesign(0.5, 7, upsampling_rate, 'sqrt');

%%
filt = conv(pulse, receivedsignal);

low = -length(freq)