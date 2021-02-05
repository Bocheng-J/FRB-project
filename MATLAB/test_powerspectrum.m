close all
clear all
clc

fs = 2000;
t = 0:1/fs:1-1/fs;
Vs1 = sin(2*pi*100*t);

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

% plot Vs1 power spectrum
powerVs1 = (abs(fftVs1).^2)/N;
figure;
plot(fax_Hz, 10*log10(abs(powerVs1)));
title('power');

figure;
pspectrum(Vs1,fs,'power');
title('power 1');

