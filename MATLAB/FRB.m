% simulation model of FRB monitoring system
clear all
close all
clc

%% Generate original signal Vs1

fs = 256e6;
t = 0:1/fs:0.000004-1/fs;
%Vs1 = chirp(t,50e6,0.000004,40e6);
Vs1 = sin(2*pi*50e6*t)+sin(2*pi*55e6*t);

N = length(Vs1);
N_2 = ceil(N/2);

fax_Hz = (0:N-1)*(fs/N);
fftVs1 = fft(Vs1);

figure;
subplot(3,1,1);
plot(t*1e6,Vs1);
title('Vs1 time-domain');
xlabel('time/us');ylabel('voltage/v');

subplot(3,1,2);
plot(fax_Hz(1:N_2)/1e6,abs(fftVs1(1:N_2)));
title('Vs1 frequency-domain')
xlabel('frequency/MHz');ylabel('magnitude');

% subplot(3,1,3);
% plot(fax_Hz(1:N_2)/1e6,abs(fftVs1(1:N_2)));
% title('Vs1 frequency-domain')
% xlabel('frequency/MHz');ylabel('magnitude');

%% ADC 
adcBitWidth = 12;
Vs1 = Vs1+abs(min(Vs1));
Vs1 = round(Vs1*(adcBitWidth-1));
fftVs1 = fft(Vs1);
% 
figure;
subplot(2,1,1);
plot(t,Vs1);
title('Vs1 time-domain 12bit quantization');
xlabel('time/us');ylabel('magnitude');

subplot(2,1,2);
plot(fax_Hz(1:N_2)/1e6,abs(fftVs1(1:N_2)));
title('Vs1 frequency-domain 12bit quantization')
xlabel('frequency/MHz');ylabel('magnitude');

%% Tunning
% Down conversion mixer
flo = 39e6;
LO = sin(2*pi*flo*t);

tunnedVs1 = Vs1.*LO;

ffttunnedVs1 = fft(tunnedVs1);

figure;
plot(fax_Hz(1:N_2)/1e6,abs(ffttunnedVs1(1:N_2)));

% lowpass filter
lowpassVs1 = lowpass(tunnedVs1,40e6,fs);

figure;
plot(t,lowpassVs1);

figure;
fftlowpassVs1 = fft(lowpassVs1);
plot(fax_Hz(1:N_2)/1e6,abs(fftlowpassVs1(1:N_2)));



