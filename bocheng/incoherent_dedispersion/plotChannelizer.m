function [] = plotChannelizer(sig,t,D,fs,flag)
% plot channelizer output in time domain figure and magnitude spectrum
% sig: channelizer output
% t: time array before channelizer
% fs: sampling frequency before channelizer
% flag: 0 for time domain figure; 1 for magnitude spectrum
figure;

 for i = 1:D                            % calculate the y axis range for magnitude spectrum
   fftsig = fft(sig(i,:));
   fft_max(1,i) = max(abs(fftsig));
   sig_max(1,i) = max(real(sig(i,:)));
   sig_min(1,i) = min(real(sig(i,:)));
   y_lim_spectrum = 1.2*max(fft_max);    
   y_max_sig = 1.2*max(sig_max);
   y_min_sig = 1.2*min(sig_min);
 end

switch flag
    case 0
       sgtitle('channelizer output (time domain)'); 
       for i=1:D
           subplot(4,D/4,i);
      % plot(downsample(t,D)*1e6,real(sig(i,:)));grid on;   % time in us
           plot(downsample(t,D),real(sig(i,:)));grid on;  % time in s
           ylim([y_min_sig y_max_sig]);
           title(['channel',num2str(i)]);
           xlabel('t/s');ylabel('voltage/v');
       end
    case 1
        sgtitle('channelizer output (magnitude spectrum)');
        for i=1:D
            N = length(sig(i,:));
            %N_2 = ceil(N/2);
            fax_Hz = (0:N-1)*(fs/D/N);
            fftsig = fft(sig(i,:));
            subplot(4,D/4,i);
            plot(fax_Hz,abs(fftsig));                          % frequency in Hz
            plot(fax_Hz/1e6,abs(fftsig));                      % frequency in MHz         
            ylim([0 y_lim_spectrum]);
            title(['channel',num2str(i)]);
            xlabel('frequency/Hz');ylabel('magnitude');
        end
end
       
end

