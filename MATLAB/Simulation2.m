close all; clear all;

%% Generate signal Vs1
fs = 7000; % sampling frequency in Hz
f = 1000; % frequency of signal Vs1 in Hz
A = 1; % amplitude of signal Vs1 in watts
t = 0:1/fs:1-1/fs;
Vs1 = A*sin(2*pi*f*t);

figure;
plot(t, Vs1);
title('Time-Domain signal');
fftVs1 = fft(Vs1)./numel(Vs1);
fftVs1 = fftshift(fftVs1);
fr = fs/2*linspace(-1,1,fs);

figure;
plot(fr, abs(fftVs1));
title('magnitude FFT of sine');
xlabel('Frequency (Hz)');
ylabel('magnitude');

%% Amplifier A1
G1 = 10;                        % voltage gain of A1 in dB
Ph1 = 45;                      % phase introduced by A1
Vn1 = 1;                        % noise introduced by A1 in dB

Ph1 = Ph1*pi/180; 
Vs1 = A*sin(2*pi*f*t+Ph1);
N = size(Vs1);
sigPower_Lin = sum(Vs1.^2)/N(2);       % signal power in watt
sigPower = 10*log10(sigPower_Lin);  % signal power in dB
noisePower_Lin = 10^(Vn1/10);           % noise power in watt
noise = sqrt(noisePower_Lin)*randn(1,N(2));    % generate noise signal
A1_Out = Vs1+noise;                 % add noise to signal
G1_Lin = 10^(G1/20); % convert voltage gain from dB to linear
A1_Out = A1_Out.* G1_Lin;           % Output of amplifier A1

figure;
plot(t,A1_Out);            % times 1000 to change the unit from s to ms


fftA1 = fft(A1_Out)./numel(A1_Out);
fftA1 = fftshift(fftA1);
figure;
plot(fr, abs(fftA1));
title('magnitude FFT of noise');
xlabel('Frequency (Hz)');
ylabel('magnitude');





