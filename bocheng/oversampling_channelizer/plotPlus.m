function [] = plotPlus(sig,t,fs,toptitle)
% plot time domain figure, magnitude spectrum and power spectrum in one
% figure
% sig: input signal
% t: x-axis of time domain figure
% fs: sampling frequency of input signal
% toptitle: top title of the figure

figure;
sgtitle(toptitle);
subplot(3,1,1);
%plot(t*1e6,real(sig)); grid on;                        % time in us
plot(t,real(sig)); grid on;                        % time in s
title('time domain');
xlabel('t/s'); ylabel('voltage/v');

subplot(3,1,2);
N = length(sig);
N_2 = ceil(N/2);
fax_Hz = (0:N-1)*(fs/N);
fftsig = fft(sig);
%plot(fax_Hz(1:N_2)/1e6,abs(fftsig(1:N_2))); grid on;           % frequency in MHz
plot(fax_Hz(1:N_2),abs(fftsig(1:N_2))); grid on;           % frequency in Hz
title('magnitude spectrum');
xlabel('frequency/Hz'); ylabel('magnitude');

subplot(3,1,3);
pow = abs(fftsig).^2;
%plot(fax_Hz(1:N_2)/1e6,10*log10(pow(1:N_2))); grid on;      % frequency in MHz
plot(fax_Hz(1:N_2),10*log10(pow(1:N_2))); grid on;      % frequency in Hz
title('power spectrum');
xlabel('frequency/Hz'); ylabel('power/dB');
end

