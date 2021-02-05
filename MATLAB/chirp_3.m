%% program initialization 
clc;
close all;
clear;

fs = 256*1e6;                                      % sample frequency
t = 0:1/fs:0.000004-1/fs;                          % time length : 4us


%% generate chirp signal
<<<<<<< HEAD
X = chirp(t,100e6,0.000002,50e6);                 % linear chirp: at time=0,frequency = 100MHz; at time=2us, frequency = 50MHz.
=======
X = chirp(t,50e6,0.000004,40e6);                 % linear chirp: at time=0,frequency = 100MHz; at time=2us, frequency = 50MHz.
>>>>>>> xiaochen
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
str_title2 = sprintf('Signal Spectrum: S1');
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
% to be finished


%% add noise and RFI
% Vn1 = 0;                                        % noise introduced by A1 in dB
% noisePower_Lin = 10^(Vn1/10);                   % noise power in watt
%noise = 0.1*randn(1,N);        % generate noise signal
noise = rand(1,N)*0.1;

f_rfi = 45e6;                         % RFI frequency: 60M hz
S_rfi = 0.1*sin(2*pi*f_rfi*t);       % RFI at 60 MHz
Vs1 = X + noise + S_rfi;

SNR = snr(X, noise+S_rfi);           % calculate SNR
disp_snr = sprintf('SNR = %f', SNR);
disp(disp_snr);

figure(3);
subplot(3,1,1);
plot(t/1e-6,Vs1);
str_title6 = sprintf('Signal S1+noise+RFI, fs = %d MHz', fs/1e6);
title(str_title6);
xlabel('t/us'); ylabel('volt/v'); 

subplot(3,1,2);
X_mags_2 = abs(fft(Vs1));
plot((fax_Hz(1:N_2))/1e6, X_mags_2(1:N_2));
str_title5 = sprintf('Signal Spectrum: S1+noise+RFI');
title(str_title5);
xlabel('Frequency/MHz'); ylabel('Magnitude');


%% Power spectrum: Vs1+noise+RFI 

subplot(3,1,3);
Pout_2 = (fft(Vs1).^2);
Pout_2_db = 10*log10(abs(Pout_2));
plot((fax_Hz(1:N_2))/1e6, Pout_2_db(1:N_2))                  % power spectrum calculated by FFT function
str_title6 = sprintf('Signal Power Spectrum: S1+noise+RFI, SNR=%f dB',SNR);
title(str_title6);
xlabel('Frequency/MHz'); ylabel('Power/dB');

%% Tuning 
% Down conversion mixer
flo = 30e6;
LO = sin(2*pi*flo*t);

tunnedVs1 = Vs1.*LO;

<<<<<<< HEAD
posVs1 = Vs1+abs(min(Vs1));                                 % move origin signal above zero    
% figure;
% plot(t,posVs1);
% title('lift signal amplitude baseband');
% xlabel('t/us'); ylabel('Frequency/MHz'); 

figure;
subplot(3,1,1);
quantizedVs1 = round(posVs1*(bitWidth-1));                  % after quantization (in time domain)
plot(t,quantizedVs1);
title('Quantization result: time domain')
xlabel('t/us'); ylabel('Frequency/MHz'); 

% r=ceil(y*(Bit_Width-1));      % 
% r=floor(y*(Bit_Width-1));     % ?
%r=round(y*(Bit_Width-1));       % 

fftquantizedVs1 = fft(quantizedVs1);
S_AftQuan = abs(fftquantizedVs1(1:N_2));
S_AftQuan(1,1) = 0;

subplot(3,1,2);
plot((fax_Hz(1:N_2))/1e6, S_AftQuan);                            % after quantization (in frequency domain)
title('Quantization result: frequency domain')
xlabel('Frequency/MHz'); ylabel('Magnitude');

subplot(3,1,3);
Pout_AftQuan = (S_AftQuan.^2);
Pout_db_AftQuan = 10*log10(abs(Pout_AftQuan));
plot((fax_Hz(1:N_2))/1e6, Pout_db_AftQuan(1:N_2))                  % power spectrum calculated by FFT function

SNR_AftQuan = 6.02*bitWidth + 1.76;                                         % calculate SNR 
%Q_level = max(posVs1)-min(posVs1);                                          % calculate quantization level
Q_level = 2 + 0.2 + 0.2;
Q_noise = 2.4/(12)^0.5;                                                 % calculate Quantization noise
Q_noise_db = 10*log10(Q_noise);                                             
disp_qnoise = sprintf('Quantization noise = %f dB', Q_noise_db);
disp(disp_qnoise);
disp_snr_AftQuan = sprintf('SQNR(ideal) = %f dB', SNR_AftQuan);
disp(disp_snr_AftQuan);

str_title7 = sprintf('Quantization result: Power Spectrum, SNR(ideal)=%f dB', SNR_AftQuan);
title(str_title7);
xlabel('Frequency/MHz'); ylabel('Power/dB');
=======
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

%%  Quantization
bitWidth = 8;

posVs1 = lowpassVs1+abs(min(lowpassVs1));                                 % move origin signal above zero    
% figure;
% plot(t,posVs1);
% title('lift signal amplitude baseband');
% xlabel('t/us'); ylabel('Frequency/MHz'); 

figure;
subplot(3,1,1);
quantizedVs1 = round(posVs1*(bitWidth-1));                  % after quantization (in time domain)
plot(t,quantizedVs1);
title('Quantization result: time domain')
xlabel('t/us'); ylabel('Frequency/MHz'); 

% r=ceil(y*(Bit_Width-1));      % 
% r=floor(y*(Bit_Width-1));     % ?
%r=round(y*(Bit_Width-1));       % 

fftquantizedVs1 = fft(quantizedVs1);
S_AftQuan = abs(fftquantizedVs1(1:N_2));
S_AftQuan(1,1) = 0;
>>>>>>> xiaochen

subplot(3,1,2);
plot((fax_Hz(1:N_2))/1e6, S_AftQuan);                            % after quantization (in frequency domain)
title('Quantization result: frequency domain')
xlabel('Frequency/MHz'); ylabel('Magnitude');

<<<<<<< HEAD
%% Tuning 
=======
subplot(3,1,3);
Pout_AftQuan = (S_AftQuan.^2);
Pout_db_AftQuan = 10*log10(abs(Pout_AftQuan));
plot((fax_Hz(1:N_2))/1e6, Pout_db_AftQuan(1:N_2))                  % power spectrum calculated by FFT function

SNR_AftQuan = 6.02*bitWidth + 1.76;                                         % calculate SNR 
%Q_level = max(posVs1)-min(posVs1);                                          % calculate quantization level
Q_level = 2 + 0.2 + 0.2;
Q_noise = 2.4/(12)^0.5;                                                 % calculate Quantization noise
Q_noise_db = 10*log10(Q_noise);                                             
disp_qnoise = sprintf('Quantization noise = %f dB', Q_noise_db);
disp(disp_qnoise);
disp_snr_AftQuan = sprintf('SQNR(ideal) = %f dB', SNR_AftQuan);
disp(disp_snr_AftQuan);

str_title7 = sprintf('Quantization result: Power Spectrum, SNR(ideal)=%f dB', SNR_AftQuan);
title(str_title7);
xlabel('Frequency/MHz'); ylabel('Power/dB');
>>>>>>> xiaochen


