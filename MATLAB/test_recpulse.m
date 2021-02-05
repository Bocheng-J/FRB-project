clear all
close all
clc

fs = 1000;
t = -0.5:1/fs:0.5-1/fs;
T = 0.2;
Vs1 = rectpuls(t, 0.2);
figure;
%plot(t,Vs1);

spectrogram(Vs1,512,500,1024,1E3,'yaxis')

figure;
plot(t, Vs1, 'k');

N = length(Vs1);                                                            % length of signal Vs1
fftVs1 = fft(Vs1);
fftVs1 = fftshift(fftVs1);
fax_Hz = (0:N-1)*(fs/N);                                                    % frequency axis
figure;
plot(fax_Hz, abs(fftVs1));                        % plot single sided frequency spectrum
title('Signal Vs1 (frequency-domain)');
xlabel('frequency/Hz'); ylabel('magnitude'); grid on;