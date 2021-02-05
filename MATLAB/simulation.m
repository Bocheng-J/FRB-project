close all; clear all;

%% Generate signal Vs1
% parameter
fs = 8000;                                                                  % sampling frequency in Hz
f = 100;                                                                   % frequency of signal Vs1 in Hz
A = 5;                                                                      % amplitude of signal Vs1 in watts
t = 0:1/fs:0.1-1/fs;                                                          % length of x-axis

% generate signal Vs1
Vs1 = A*sin(2*pi*f*t);                                                      

% plot Vs1 in time domain
figure;
subplot(2,1,1);
plot(t, Vs1);
title('Signal Vs1 (time-domain)');
xlabel('t/s'); ylabel('volt/v'); grid on;

% plot Vs1 in frequency domain
N = length(Vs1);                                                            % length of signal Vs1
fftVs1 = fft(Vs1);                                                         
fax_Hz = (0:N-1)*(fs/N);                                                    % frequency axis
subplot(2,1,2);
plot(fax_Hz(1:ceil(N/2)), abs(fftVs1(1:ceil(N/2))));                        % plot single sided frequency spectrum
title('Signal Vs1 (frequency-domain)');
xlabel('frequency/Hz'); ylabel('magnitude'); grid on;

%% Amplifier A1
% parameter
G1 = 10;                                                                    % gain in dB
Ph1 = 45;                                                                   % added phase in 
Vn1 = 10;                                                                  % added noise in dB

% add phase
Ph1 = Ph1*pi/180;                                                           % convert degree to rad
%Vs1 = A*sin(2*pi*f*t+Ph1);                                                  % add phase to signal Vs1
Vs1 = circshift(Vs1, 45);
% add noise
sigPower = sum(Vs1.^2)/N;                                                   % calculate signal Vs1 power in watt
sigPower = 10*log10(sigPower);                                              % convert signal power from watt to dB
noisePower = 10^(Vn1/10);                                                   % calculate noise power in watt
ampNoise = sqrt(noisePower)*randn(1,N);                                     % generate noise signal
SNR = snr(Vs1, ampNoise);                                                   % calculate SNR in dB
Amp1 = Vs1+ampNoise;                                                        % add amplifier noise to signal Vs1

% add gain
G1 = 10^(G1/20);                                                            % convert voltage gain from dB to linear
Amp1 = Amp1.* G1;                                                           % amplify signal Vs1

% plot output of amplifier A1 in time domain
figure;
subplot(3,1,1);
plot(t,Amp1); hold on;
% plot(t,ampNoise); hold on;
% plot(t,Vs1); hold on;
title('Output of amplifier A1(time-domain)');
xlabel('t/s'); ylabel('volt/v'); grid on;

% plot output of amplifier A1 in frequency domain
fftAmp1 = fft(Amp1);                                                         
subplot(3,1,2);
plot(fax_Hz(1:ceil(N/2)), abs(fftAmp1(1:ceil(N/2))));                         % plot single sided frequency spectrum
title('Output of amplifier A1(frequency-domain)');
xlabel('frequency/Hz'); ylabel('magnitude'); grid on;

% plot power spectrum of output of amplifier A1
absfftA1 = abs(fftAmp1);                                                      
powerAmp1 = absfftA1.^2/(N);
subplot(3,1,3);
plot(fax_Hz(1:ceil(N/2)) ,10*log10(powerAmp1(1:ceil(N/2)))); grid on;       % plot power spectrum in dB
title('Output of amplifier A1(power spectrum)');
xlabel('frequency/Hz'); ylabel('power/dB'); grid on;

%% Add dispersion
figure;
Amp1 = ifft(fftAmp1);
plot(Amp1);



