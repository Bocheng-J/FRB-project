clear all
close all
clc

%%
fs = 1000;
t = (0:1/fs:1-1/fs)';

Vs1 = chirp(t,220,t(end),0);
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




