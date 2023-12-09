%% Setup
clear; clc; close all;

L = 32; %processing gain/spreading gain

sync_len = 20; % in samples
fs = 200; %MHz
upsampling_rate = 12;
new_preamble = 0;

if(new_preamble)
    disp("New Preamble")

    freq_sync = ones(1, sync_len); %2 us of frame sync for the receiver to get PLL lock
    t_sync = floor(2.*rand(1, sync_len)); %Pseudo-random
    frame_sync = floor(2.*rand(1, sync_len)); %Pseudo-random
    gc1 = generate_gold_code(32);
    save("gold_codes.mat", "gc1");
    save("frame_sync", "frame_sync")
    save("t_sync", "t_sync")
    save("freq_sync", "freq_sync")
else
    load gold_codes
    load freq_sync;
    load t_sync.mat;
    load frame_sync.mat;
end

preamble = horzcat(freq_sync, t_sync, frame_sync);

if(mod(length(frame_sync), 2)) error("NOT MOD 2"); end

%% Image to QAM

image = imread('images/shannon1440.bmp');
dims = size(image);
im = double(reshape(image, 1, dims(1)*dims(2)));

user1 = im(1:100);
user2 = im(101:200);
n_frames = 10;

im = reshape(user1, length(user1)/n_frames, n_frames);
bits1 = horzcat(preamble, im(:,1)', frame_sync, im(:,2)', frame_sync, im(:,3)', ...
    frame_sync, im(:,4)', frame_sync, im(:,5)', frame_sync, im(:,6)', frame_sync, ...
    im(:,7)', frame_sync, im(:,8)', frame_sync, im(:,9)', frame_sync, im(:,10)');

im = reshape(user2, length(user2)/n_frames, n_frames);
bits2 = horzcat(preamble, im(:,1)', frame_sync, im(:,2)', frame_sync, im(:,3)', ...
    frame_sync, im(:,4)', frame_sync, im(:,5)', frame_sync, im(:,6)', frame_sync, ...
    im(:,7)', frame_sync, im(:,8)', frame_sync, im(:,9)', frame_sync, im(:,10)');

bits1 = string(reshape(bits1, 2, length(bits1)/2));
bits2 = string(reshape(bits2, 2, length(bits2)/2));

qam_bits1 = strings(1, length(bits1));
qam_bits2 = strings(1, length(bits2));
for i = 1:length(bits1)
    qam_bits1(i) = strjoin(bits1(:,i),'');
    qam_bits2(i) = strjoin(bits2(:,i),'');
end

qam_bits1 = bin2dec((qam_bits1))+1;
qam_bits2 = bin2dec((qam_bits2))+1;

d = 2;
options = [0.5+0.5j, -0.5-0.5j
           0.5-0.5j, -0.5+0.5j].*d;

qam_points1 = options(qam_bits1);
qam_points2 = options(qam_bits2);

qam_points1 = repmat(qam_points1, L, 1);
qam_points2 = repmat(qam_points2, L, 1);

qam_points1 = reshape(qam_points1, 1, L*length(bits1));
qam_points2 = reshape(qam_points2, 1, L*length(bits2));

% Assign unique Gold Codes for each user
gc2 = gc1(2, :);
gc1 = gc1(1, :);

%Spread signal based on assigned Gold Code
spread1 = spread_signal_gold(qam_points1, gc1, 1);
spread2 = spread_signal_gold(qam_points2, gc2, 1);

summed = spread1 + spread2; 

qam_preamble = summed(1:length(preamble)/2*L);
qam_frame_sync = qam_preamble(end-length(frame_sync)*L/2+1:end);
save('qam_pre', 'qam_preamble');
save('qam_fsync', 'qam_frame_sync');

%% Transmission

bits_up = upsample(summed, upsampling_rate);
pulse = rcosdesign(0.3, 40, upsampling_rate, 'sqrt');
transmitsignal = conv(bits_up, pulse);
transmitsignal = transmitsignal./max(abs(transmitsignal));
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

figure;
subplot(3, 1, 1); hold on
plot(real(spread1(1:100)));
plot(imag(spread1(1:100))); title("CH1"); hold off
subplot(3, 1, 2); hold on
plot(real(spread2(1:100)));
plot(imag(spread2(1:100)));title("CH2"); hold off
subplot(3, 1, 3); hold on
plot(real(summed(1:100)));
plot(imag(summed(1:100)));title("Combined"); hold off


%% Gold Code Generation

function gold_codes = generate_gold_code(N)
    user1 = randsrc(1,5,[0 1]); user2 = randsrc(1,5,[0 1]);
    m_seq1 = zeros(1, N); m_seq2 = zeros(1, N);
    for j = 1:N
        m_seq1(j) = user1(1);
        if(user1(1) == user1(4))
            temp3 = 0;
        else
            temp3 = 1;
        end
        user1(1) = user1(2);
        user1(2) = user1(3);
        user1(3) = user1(4);
        user1(4) = user1(5);
        user1(5) = temp3;
    end
    for i = 1:N
        m_seq2(i) = user2(1);
        if(user2(1) == user2(2))
            temp3 = 0;
        else
            temp3 = 1;
        end
        if(user2(4) == temp3)
            temp2 = 0;
        else
            temp2 = 1;
        end
        if(user2(5) == temp2)
            temp3 = 0;
        else
            temp3 = 1;
        end
        user2(1) = user2(2);
        user2(2) = user2(3);
        user2(3) = user2(4);
        user2(4) = user2(5);
        user2(5) = temp3; 
    end
    gold_codes = zeros(2, N);
    for user = 1:2
        temp = m_seq2(1);
        
        %Shift by 1 bit
        for i = 1:N-1
            m_seq2(i) = m_seq2(i+1);
        end
        m_seq2(N) = temp;
        
        %Code Generation
        codes = zeros(1, N);
        for i = 1:N
            codes(i) = xor(m_seq1(i), m_seq2(i));
        end
        
        gold_codes(user, :) = codes;
    end
    
    %Normalise
    for r = 1:2
        for c = 1:N
            if gold_codes(r,c) == 0
                gold_codes(r,c) = -1;
            end
        end
    end
end

function signal_spread = spread_signal_gold(data, gold_code, spreading_gain)
    repeated_gold_code = repmat(gold_code, 1, ceil(length(data) / length(gold_code)));
    signal_spread = data .* repeated_gold_code(1:end) * spreading_gain;
end
