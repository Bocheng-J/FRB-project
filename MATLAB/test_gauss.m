% Plot a 50 kHz Gaussian RF pulse with 60% bandwidth, sampled at a rate of
% 10 MHz. Truncate the pulse where the envelope falls 40 dB below the peak.

fs = 256e6;
tc = gauspuls('cutoff',50e6,0.6,[],-40); 
t = -tc/2:1/fs:tc/2;  
yi = gauspuls(t,60e6,0.6); 

plot(t,yi);


% tc = gauspuls('cutoff',5e6,0.5,[],-40);
% t = (-tc) : 1e-8 : (tc);
% yi = gauspuls(t,5e6,0.5);
% %Addition of the zeros, before I implement the shift
% yi1 = [yi; zeros(1,100)];
% %This is where the error occurs



% %% generate Vs1
% A = 1;
% fs = 8000;
% f1 = 250;
% f2 = 780;
% f3 = 1700;
% t = 0:1/fs:0.1-1/fs; 
% Vs1 = A*(5*sin(2*pi*f1*t)+10*sin(2*pi*f2*t)+15*sin(2*pi*f3*t));
% 
% figure;
% subplot(2,1,1);
% plot(t, Vs1);
% title('Signal Vs1 (time-domain)');
% xlabel('t/s'); ylabel('volt/v'); grid on;
% 
% % plot Vs1 in frequency domain
% N = length(Vs1);                                                            % length of signal Vs1
% fftVs1 = fft(Vs1);                                                         
% fax_Hz = (0:N-1)*(fs/N);                                                    % frequency axis
% subplot(2,1,2);
% plot(fax_Hz(1:ceil(N/2)), abs(fftVs1(1:ceil(N/2))));                        % plot single sided frequency spectrum
% title('Signal Vs1 (frequency-domain)');
% xlabel('frequency/Hz'); ylabel('magnitude'); grid on;
% 
% %% add dispersion
% 
% % divide the frequency band into a number of frequency bins
% yfilt = filter1('bp',Vs1,'fc',[0 20],'fs',Fs);
% 
% 
% 
% 
% 
% 
% 
% % bin0 = bandpass(Vs1,[1 2000],fs);
% % bin1 = bandpass(Vs1,[501 1000],fs);
% % bin2 = bandpass(Vs1,[1001 1500],fs);
% % bin3 = bandpass(Vs1,[1501 2000],fs);



