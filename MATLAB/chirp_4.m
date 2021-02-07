%% program initialization 
clc;
close all;
clear;

fs = 256*1e6;                                      % sample frequency
t = 0:1/fs:0.000004-1/fs;                          % time length : 4us


%% generate chirp signal
%X = chirp(t,50e6,0.000004,40e6);                 % linear chirp: at time=0,frequency = 100MHz; at time=2us, frequency = 50MHz.
A = 1;
X = [zeros(1,1024), A*chirp(t,50e6,0.000004,40e6)];
refX = A*chirp(t,50e6,0.000004,40e6);
t1 = 0:1/fs:0.000004*2-1/fs;
figure(1);
subplot(2,1,1);
plot(t1/1e-6,X);
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
% [S,F,T,P]= spectrogram(X,window,noverlap,nfft,fs,'yaxis');
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
noise = randn(1,N)*0.1;

f_rfi = 45e6;                         % RFI frequency: 60M hz
S_rfi = 0.1*sin(2*pi*f_rfi*t1);       % RFI at 60 MHz
Vs1 = X + noise + S_rfi;

SNR = snr(X, noise+S_rfi);           % calculate SNR
disp_snr = sprintf('SNR = %f', SNR);
disp(disp_snr);

figure(3);
subplot(3,1,1);
plot(t1/1e-6,Vs1);
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
flo = 30e6;                                                      % tuning frequency
LO = sin(2*pi*flo*t1);                                           % sine wave used in tuning
refLO = sin(2*pi*flo*t);                                         % reference template signal for matched filter

tunnedVs1 = Vs1.*LO;                                             % tuning
ffttunnedVs1 = fft(tunnedVs1);    

tunnedRefVs1 = refX.*refLO;                                      % tune the reference template signal

figure;
subplot(3,1,1);
plot(fax_Hz(1:N_2)/1e6,abs(ffttunnedVs1(1:N_2)));
str_title8 = sprintf('Tune to 20 MHz');
title(str_title8);
xlabel('Frequency/MHz'); ylabel('Magnitude');

lowpassVs1 = lowpass(tunnedVs1,40e6,fs);                       % lowpass filter
lowpassRefVs1 = lowpass(tunnedRefVs1,40e6,fs);

subplot(3,1,3);
plot(t1/1e-6,lowpassVs1);
title('Low-pass result in time domain');
xlabel('t/us'); ylabel('volt/v'); 

subplot(3,1,2);
fftlowpassVs1 = fft(lowpassVs1);
plot(fax_Hz(1:N_2)/1e6,abs(fftlowpassVs1(1:N_2)));
title('Low-pass spectrum');
xlabel('Frequency/MHz'); ylabel('Magnitude');

fftlowpassRefVs1 = fft(lowpassRefVs1);

% figure;
% plot(t,lowpassRefVs1);
% title('lowpass RefVs1 in time domain');
% 
% figure;
% plot(fax_Hz(1:N_2)/1e6,abs(fftlowpassVs1(1:N_2)));
% title('lowpass RefVs1 spectrum');


%% matched filter
Ar = 1;
refVs1 = lowpassRefVs1;                                                 % matched filter template
energy_refVs1 = sum(refVs1.^2);
refVs1 = refVs1/sqrt(energy_refVs1);                                    % normalization

figure;
subplot(3,1,1);
plot(fax_Hz(1:N_2)/1e6,abs(fftlowpassVs1(1:N_2)));
title('Matched filter template');
xlabel('Frequency/MHz'); ylabel('Magnitude');

corr = conv(lowpassVs1,fliplr(refVs1));                                % matched filtering: correlation
subplot(3,1,2);
plot(corr);
title('Correlation result');

[tmp, Tmax] = max(corr);                                               % find the most matched point
Tx_hat = Tmax-length(refVs1)                                           % calculate time delay of the most matched point

N = length(corr);
N_2 = ceil(N/2);
fax_Hz = (0:N-1)*(fs/N);
fftcorr = fft(corr);
subplot(3,1,3);
plot(fax_Hz(1:N_2)/1e6,abs(fftcorr(1:N_2)));                         
title('Spectrum after matched filter');
xlabel('Frequency/MHz'); ylabel('Magnitude');


%%  Quantization
bitWidth = 8;

posVs1 = corr(length(corr)-Tx_hat+1:end)+abs(min(corr(length(corr)-Tx_hat+1:end)));                                 % move origin signal above zero    
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
Nq = length(quantizedVs1);
Nq_2 = ceil(Nq/2);
S_AftQuan = abs(fftquantizedVs1(1:Nq_2));
S_AftQuan(1,1) = 0;

faxq_Hz = (0:Nq-1)*fs/Nq;


subplot(3,1,2);
plot((faxq_Hz(1:Nq_2))/1e6, S_AftQuan);                            % after quantization (in frequency domain)
title('Quantization result: frequency domain')
xlabel('Frequency/MHz'); ylabel('Magnitude');

subplot(3,1,3);
Pout_AftQuan = (S_AftQuan.^2);
Pout_db_AftQuan = 10*log10(abs(Pout_AftQuan));
plot((faxq_Hz(1:Nq_2))/1e6, Pout_db_AftQuan(1:Nq_2))                  % power spectrum calculated by FFT function

SNR_AftQuan = 6.02*bitWidth + 1.76;                                         % calculate SNR 
%Q_level = max(posVs1)-min(posVs1);                                          % calculate quantization level
% Q_level = 2 + 0.2 + 0.2;
% Q_noise = 2.4/(12)^0.5;                                                 % calculate Quantization noise
% Q_noise_db = 10*log10(Q_noise);                                             
% disp_qnoise = sprintf('Quantization noise = %f dB', Q_noise_db);
% disp(disp_qnoise);
disp_snr_AftQuan = sprintf('SQNR(ideal) = %f dB', SNR_AftQuan);
disp(disp_snr_AftQuan);

str_title7 = sprintf('Quantization result: Power Spectrum, SNR(ideal)=%f dB', SNR_AftQuan);
title(str_title7);
xlabel('Frequency/MHz'); ylabel('Power/dB');

% 0207-1823

