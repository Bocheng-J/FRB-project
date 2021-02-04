%% program initialization 
clc;
close all;
clear;

fs = 256*1e6;                                      % sample frequency
t = 0:1/fs:0.000004-1/fs;                          % time length : 4us

%% generate chirp signal
 X = chirp(t,100e6,0.000002,50e6);                 % linear chirp: at time=0,frequency = 100MHz; at time=2us, frequency = 50MHz.
figure(1);
subplot(2,1,1);
plot(t/1e-6,X);
str_title1 = sprintf('Signal S1, fs = %d MHz', fs/1e6);
title(str_title1);
xlabel('t/us'); ylabel('volt/v'); 

%% FFT
figure(2);
 subplot(2,1,1);
 N= length(X);
fax_bins = [0 : N-1];
fax_Hz = fax_bins*(fs/N);
X_mags = abs(fft(X));
N_2 = ceil(N/2);
plot((fax_Hz(1:N_2))/1e6, X_mags(1:N_2));
str_title2 = sprintf('Signal S1 Spectrum');
title(str_title2);
xlabel('Frequency/MHz'); ylabel('Magnitude');

%% Time vs Frequency plot
figure(1);
subplot(2,1,2);
window = hamming(128);  			% haming window
noverlap = 120; 					%overlap length
nfft = 2^nextpow2(length(window)); 	% dft number
spectrogram(X,window,noverlap,nfft,fs,'yaxis');
str_title3 = sprintf('Signal S1: Frequency vs Time');
title(str_title3);
xlabel('t/us'); ylabel('Frequency/MHz'); 

%% power spectrum
figure(2);
subplot(2,1,2);
Pout = (fft(X).^2);
Pout_db = 10*log10(abs(Pout));
plot((fax_Hz(1:N_2))/1e6, Pout_db(1:N_2))                  %% power spectrum calculated by FFT function
str_title4 = sprintf('Signal S1 Power Spectrum');
title(str_title4);
xlabel('Frequency/MHz'); ylabel('Power/dB');

%% channelization








